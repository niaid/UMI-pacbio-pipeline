"""
This module executes removes fake UMIs using a network process based on genetic distance of inserts
and a maximum edit distance between UMIs of 1.
"""
import sys
import time
import os
import subprocess
import itertools
import copy
from joblib import Parallel, delayed
import networkx as nx
from Bio import SeqIO
from tqdm import tqdm


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
    header_string = header_string.split("_")[-2:]
    header_string = "_".join(header_string)
    return header_string


def hamming_distance(string1, string2):
    """
    The computes hamming distance of two UMIs

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
            return False
    return True


def check_gen_dist(sample_tuple, fasta_reads):
    """
    Checks if a reads in a tuple of gen dist 0

    Parameters:

    sample_tuple: tuple
        tuples in a read to check

    fasta_reads: dict
        keys are read number and values are the read

    Return:

    sample_tuple or None:
        returns sample tuple if the gen distance is 0.
        Else it returns None
    """

    i = sample_tuple[0]
    j = sample_tuple[1]

    if hamming_distance(str(fasta_reads[i]), str(fasta_reads[j])):
        return sample_tuple


def build_gen_dist(fasta_reads, n_jobs):
    """
    Computes genetic distance between sequences. Keeps reads of gen dist 0

    Parameters
    ----------
    fasta_reads : dict
        dict of reads from fasta file

    n_jobs : int
        Number of jobs to run in parallel
    Returns
    -------
    gen_graph : networkx graph
        graph of genetic distance 0.

    """

    print("Calculating Hamming distance for sequences")
    tic = time.perf_counter()
    n = len(fasta_reads.keys())
    read_tuple_list = itertools.combinations(list(fasta_reads.keys()), 2)
    keeps = Parallel(n_jobs=n_jobs, batch_size=256)(delayed(check_gen_dist)(sample_tuple, fasta_reads) for sample_tuple
                                                    in tqdm(read_tuple_list, total=n*(n-1)/2))
    keeps = [x for x in keeps if x]  # removes the None ones
    print(f"Done calculating Hamming distance {time.perf_counter() - tic}")
    gen_graph = nx.Graph()
    gen_graph.add_nodes_from(list(range(len(fasta_reads))))
    gen_graph.add_edges_from(keeps)
    return gen_graph


def make_edge(sample_tuple, umi_count_dict):
    """
    function to pass to parallel processing to check if there is an edge

    Parameters:

    sample_tuple: tuple
        tuples of samples to check umi edit distance and criteria
        they are integers use to look up count_umi_dict

    umi_count_dict: dict
        key is read number
        value is tuple (UMI, count)
    Returns:

    sample_tuple or None:
        sample tuple if it satisfies the criteria, None otherwise
    """

    k = 0
    first_read = sample_tuple[0]
    second_read = sample_tuple[1]
    edit_dist = hamming_distance(umi_count_dict[first_read][0], umi_count_dict[second_read][0])
    if edit_dist == 1:
        count1, count2 = int(umi_count_dict[first_read][1]), umi_count_dict[second_read][1]
        criteria = max(int(count1), int(count2)) >= (2 ** (1 + 100 * (k + 0))) * min(int(count1), int(count2)) - 1
        if criteria:
            return sample_tuple


def make_umi_graph(component, umi_count_dict, n_jobs):
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
        number of jobs to run in parallel

    Returns
    -------
    umi_graph : networkx graph
        graph of umi edit distance 1

    """
    component = list(component)
    umi_tuple_list = itertools.combinations(list(component), 2)

    umi_edges = Parallel(n_jobs=n_jobs)(delayed(make_edge)(sample_tuple, umi_count_dict)
                                        for sample_tuple in umi_tuple_list)
    umi_edges = [x for x in umi_edges if x]  # Removes all None objects from the
    umi_set = [x[0] for x in umi_edges] + [x[1] for x in umi_edges]
    umi_set = list(set(umi_set))
    umi_graph = nx.Graph()
    umi_graph.add_nodes_from(umi_set)
    umi_graph.add_edges_from(umi_edges)
    return umi_graph


def make_merge_list(umi_graph, umi_count_dict, merge_dict):
    """
    Adds to merge dict

    Parameters
    ----------
    umi_graph : networkx graph
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

    umi_graph_copy = copy.deepcopy(umi_graph)
    while len(umi_graph_copy.nodes) > 1:
        umi_node_counts = [(x, umi_count_dict[x][1]) for x in umi_count_dict if x in set(umi_graph_copy.nodes)]
        umi_node_counts.sort(key=lambda x: int(x[1]), reverse=True)
        largest_node = umi_node_counts[0][0]
        neighbors = list(umi_graph_copy.neighbors(largest_node))
        merge_dict[largest_node] = neighbors
        umi_graph_copy.remove_nodes_from([largest_node]+neighbors)
    return merge_dict


def make_header_dict(merge_dict, fasta_file_dict, umi_count_dict):
    """
    makes header dict

    Parameters
    ----------
    merge_dict : dict
        tracks which merges to make
        key: umi with largest count
        value: list of reads to be merged
    fasta_file_dict : dict
        key is read number
        value is read in SeqIO format
    umi_count_dict : dict
        key is read number
        value is (umi, count)
    Returns
    -------
    header_dict : dict
        Which headers need to change
    umi_dict : dict
        key is old umi and count
        value is new umi and count

    remove_list : list
        list of all umi_headers that are removed
    """

    umi_dict = {}
    header_dict = {}
    remove_list = []

    print("Down-selecting UMIs")

    for merge_key in merge_dict:
        old_header = fasta_file_dict[merge_key].id
        new_count = sum([int(umi_count_dict[x][1]) for x in [merge_key] + merge_dict[merge_key]])
        new_header = "_".join(old_header.split('_')[:-1]) + '_' + str(new_count)
        for item in [merge_key] + merge_dict[merge_key]:
            umi_dict[strip_umi(fasta_file_dict[item].id)] = "_".join(new_header.split('_')[-2:])
        for item in merge_dict[merge_key]:
            remove_list.append("_".join(umi_count_dict[item]))
        header_dict[merge_key] = new_header
    return header_dict, umi_dict, remove_list


def main(project, file_path, n_jobs):
    """
    See module doc string
    """
    fasta_path = project + "sgs-prelim/mafft/" + file_path + ".fasta"
    fasta_file = list(SeqIO.parse(fasta_path, "fasta"))
    fasta_reads = [read.seq for read in fasta_file]
    fasta_reads = dict(zip(range(len(fasta_file)), fasta_reads))
    fasta_file_dict = dict(zip(range(len(fasta_file)), fasta_file))
    umi_count_dict = {x: (strip_umi(fasta_file_dict[x].id).split('_')[0],
                          strip_umi(fasta_file_dict[x].id).split('_')[1])
                      for x in fasta_file_dict}

    gen_graph = build_gen_dist(fasta_reads, n_jobs)

    merge_dict = {}
    for component in nx.connected_components(gen_graph):
        if len(component) > 1:
            umi_graph = make_umi_graph(component, umi_count_dict, n_jobs)
            merge_dict = make_merge_list(umi_graph, umi_count_dict, merge_dict)

    header_dict, umi_dict, remove_nw_list = make_header_dict(merge_dict, fasta_file_dict, umi_count_dict)
    remove_list = [item for sublist in merge_dict.values() for item in sublist]

    for key in remove_list:
        del fasta_file_dict[key]
    for key in header_dict:
        fasta_file_dict[key].id = header_dict[key]
        fasta_file_dict[key].name = header_dict[key]
        fasta_file_dict[key].description = header_dict[key]

    write_folder = project + "fake-umi-curation/" + file_path
    if not os.path.exists(write_folder):
        subprocess.check_call(['mkdir', write_folder])

    print("Writing Selected UMIs")
    with open(write_folder + "/" + file_path + '_post_nw.fasta', 'w') as outfile:
        SeqIO.write(list(fasta_file_dict.values()), outfile, 'fasta')

    with open(write_folder + "/" + file_path + '_merged_list.csv', 'w') as handle:
        print("old UMI count, new UMIs", file=handle)
        for key in umi_dict:
            print(f"{key},{umi_dict[key]}", file=handle)
    
    with open(write_folder + "/" + file_path + '_removed_nw.txt', 'w') as handle:
        for removed_umi in remove_nw_list:
            print(removed_umi, file=handle)


if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    N_JOBS = int(sys.argv[3])
    main(PROJECT, FILE_PATH, N_JOBS)
