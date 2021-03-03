#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module executes fake umi removal without considering genetic distance
It only looks at edit distance one pairs.
"""

import sys
import os
import subprocess
import itertools
from joblib import Parallel, delayed
import networkx as nx
from Bio import SeqIO
from tqdm import tqdm

#pylint: disable=E1101, C0200
def strip_umi(header_string):
    """
    Parameters
    ----------
    header_string : str
        This is the header of the fasta file.

    Returns
    -------
    header_string : str
        umi_number of reads in that umi bin after 99% consensus clustering.

    """
    #header_string = header_string.split(".")[1]
    header_string = header_string.split("_")[3:5]
    header_string = "_".join(header_string)
    return header_string

def hamming_distance(string1, string2):
    """
    The computes hamming distance of two umis

    Parameters
    ----------
    string1 : str
        umi1
    string2 : str
        umi2

    Returns
    -------
    distance : int
        hamming distance

    """
    distance = 0
    for i, _ in enumerate(string1):
        if string1[i] != string2[i]:
            distance += 1
    return distance

def make_edge(sample_tuple, umi_count_dict):
    """
    function to pass to parallel processing to check if there is an edge

    Parameters:

    sample_tuple: tuple
        tuples of samples to check umi edit distance and criteria
        they are integers use to look up count_umi_dict

    Returns:

    sample_tuple or None:
        sample tuple if it satisties the criteria, None otherwise
    """
    i = sample_tuple[0]
    j = sample_tuple[1]
    edit_dist = hamming_distance(umi_count_dict[i][0], umi_count_dict[j][0])
    criteria = max(int(umi_count_dict[i][1]), int(umi_count_dict[j][1])) >= \
                2*min(int(umi_count_dict[i][1]), int(umi_count_dict[j][1])) - 1
    if edit_dist == 1 and criteria:
        return sample_tuple

def make_umi_graph(umi_count_dict, n_jobs):
    """
    Creates the umi edit distance 1 graph

    Parameters
    ----------
    component : set
        set of nodes in the connected component
    umi_count_dict : dict
        key is read number
        value is (umi, count)

    n_jobs : int
        The number of parallel jobs to process

    Returns
    -------
    umi_graph : networkx graph
        graph of umi edit distance 1

    """
    umi_graph = nx.Graph()
    umi_graph.add_nodes_from(list(umi_count_dict.keys()))
    umi_tuple_list = itertools.combinations(list(umi_count_dict.keys()), 2)

    umi_edges = Parallel(n_jobs=n_jobs)(delayed(make_edge)(sample_tuple, umi_count_dict) for sample_tuple in tqdm(umi_tuple_list))
    umi_edges = [x for x in umi_edges if x] # Removes all None objects from the
    umi_graph.add_edges_from(umi_edges)
    return umi_graph

def make_merge_list(subgraph, umi_count_dict, merge_dict):
    """
    Adds to merge dict

    Parameters
    ----------
    subgraph : networkx graph
        graph of umi edit distance 1
    umi_count_dict : dict
        key is read number
        value is (umi, count)
    merge_dict : dict
        tracks which merges to make
        key: umi with largest count
        value: list of reads to be merged

    Returns
    -------
    None.

    """

    subgraph_copy = nx.Graph(subgraph)
    while len(subgraph_copy.nodes) > 1:
        umi_node_counts = [(x, umi_count_dict[x][1])\
                           for x in umi_count_dict if x in set(subgraph_copy.nodes)]
        umi_node_counts.sort(key=lambda x: int(x[1]), reverse=True)
        largest_node = umi_node_counts[0][0]
        neighbors = list(subgraph_copy.neighbors(largest_node))
        merge_dict[largest_node] = neighbors
        subgraph_copy.remove_nodes_from([largest_node]+neighbors)
    return merge_dict

def main(project, file_path, n_jobs):
    """
    See module doc string

    Parameters
    ----------
    fasta_path : str
        path to fasta file

    Returns
    -------
    None.

    """

    fasta_path = project + "final_ccs_reads/" + file_path + "read.fasta"
    fasta_file = list(SeqIO.parse(fasta_path, "fasta"))
    fasta_file_dict = dict(zip(range(len(fasta_file)), fasta_file))
    umi_count_dict = {x:(strip_umi(fasta_file_dict[x].id).split('_')[0],\
                         strip_umi(fasta_file_dict[x].id).split('_')[1])\
                      for x in fasta_file_dict}
    print('making distance matrix')
    umi_graph = make_umi_graph(umi_count_dict, n_jobs)
    merge_dict = {}
    isolated_list = []
    for component in nx.connected_components(umi_graph):
        if len(component) > 1:
            subgraph = umi_graph.subgraph(list(component))
            merge_dict = make_merge_list(subgraph, umi_count_dict, merge_dict)
        else:
            isolated_list.append(list(component)[0])


    remove_list = [fasta_file_dict[item] for sublist in merge_dict.values() for item in sublist]
    keep_list = [fasta_file_dict[item] for item in list(merge_dict.keys()) + isolated_list]

    write_folder = project + "fake-umi-curation-nogen/" + file_path
    if not os.path.exists(write_folder):
        subprocess.check_call(['mkdir', write_folder])

    print("Writing Selected UMIs")
    with open(write_folder + "/" + file_path +'_keep_nw.fasta', 'w') as outfile:
        SeqIO.write(keep_list, outfile, 'fasta')

    with open(write_folder + "/" + file_path + '_remove_nw.fasta', "w") as handle:
        SeqIO.write(remove_list, handle, "fasta")

    with open(write_folder + "/" + file_path + "edit_dist_1_pairs.csv", "w") as handle:
        for item in umi_graph.edges:
            print(f"{fasta_file_dict[item[0]].id},{fasta_file_dict[item[1]].id}",\
                  file=handle)
if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    N_JOBS = int(sys.argv[3])
    main(PROJECT, FILE_PATH, N_JOBS)
