#!/usr/bin/env python3

import os
import sys
from collections import Counter
from itertools import combinations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import networkx as nx
import numpy as np


def get_counts(seq_dict):
    """
    Get counts of each nucleotide at each position (accepts gaps).

    Arguments:
        seq_dict (dict): Dictionary of sequences, should all be aligned and have equivalent lengths.
    Returns:
        List of Counter objects with length equal to the length of sequences.
    """
    return [Counter(''.join((seq_dict[seq][j] for seq in seq_dict))) for j in range(len(list(seq_dict.values())[0]))]


def determine_consensus(counts):
    cons = dict()
    for i, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c > 0.75:
            cons[i] = counter.most_common(1)[0][0]
    return cons


def clean_pb_errors(seq_dict, counts, cons):
    for index, counter in enumerate(counts):
        for char in counter:
            if counter[char] == 1 or counter[char] < 0.01*sum(counter.values()):
                # print(char, index)
                for seq in seq_dict:
                    if seq_dict[seq][index] == char and index in cons:
                        seq_dict[seq][index] = cons[index]
    return seq_dict


def get_largest_vars(seq_dict):
    counts = get_counts(seq_dict)
    putative_vars = []
    for index, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c < 0.90:
            putative_vars.append((index, counter.most_common(2)[1][0], 1-c, counter.most_common(2)[0][0]))
    putative_vars_sorted = list(sorted(putative_vars, key=lambda v: v[2]))
    if putative_vars_sorted:
        max_counts = max(c[2] for c in putative_vars_sorted)
        # print(max_counts)
    else:
        return []
    vs = list(filter(lambda c: c[2] > 0.9*max_counts, putative_vars_sorted))
    return vs


def get_seqs_w_fp(seq_dict, ls):
    s = dict()
    sn = dict()
    for seq in seq_dict:
        match = False
        for i, site in enumerate(ls):
            if seq_dict[seq][site[0]] == site[1]:
                match = True
                break
        if match:
            s[seq] = seq_dict[seq]
        else:
            sn[seq] = seq_dict[seq]
    return s, sn


def assess_linkage(sites, seq_dict):
    linkage_map = {seq: [0] * len(sites) for seq in seq_dict}
    for seq in seq_dict:
        for i, site in enumerate(sites):
            if seq_dict[seq][site[0]] == site[1]:
                linkage_map[seq][i] = 1
    linkage_map = np.array([linkage_map[seq] for seq in linkage_map])
    # print(linkage_map)
    # input()
    lsites_indices = assess_linkage_map(linkage_map)

    return [[sites[i] for i in cc] for cc in lsites_indices]


def compute_link(a, b):
    return sum(ai == bi for ai, bi in zip(a, b))/len(a)


def assess_linkage_map(lm):
    ns = lm.shape[1]
    nodes = list(range(ns))

    G = nx.Graph()
    G.add_nodes_from(nodes)

    for combo in combinations(list(range(ns)), 2):
        # print(compute_link(lm[:, combo[0]], lm[:, combo[1]]))
        # input()
        if compute_link(lm[:, combo[0]], lm[:, combo[1]]) > 0.9:
            G.add_edge(combo[0], combo[1])

    link_indices = list(nx.connected_components(G))
    return link_indices


def get_pcr_fp(seq_dict):
    if len(seq_dict) < 5:
        return []
    fp = []
    sites = get_largest_vars(seq_dict)
    # print(sites)
    # print(sites)
    lsites = assess_linkage(sites, seq_dict)
    # print(lsites)
    fp.extend(lsites)
    print(lsites)
    for ls in lsites:
        lseqs, lseqs_n = get_seqs_w_fp(seq_dict.copy(), ls)
        if len(lseqs) == len(seq_dict):
            break
        fp.extend(get_pcr_fp(lseqs))
        fp.extend(get_pcr_fp(lseqs_n))
    return fp


def flatten_fingerprint(fp_list, par):
    l = []
    pari = []
    for fp in fp_list:
        if isinstance(fp, tuple):
            # print(fp)
            pari.append((fp[0], fp[1]))
            continue
        if isinstance(fp, list):
            pari = par + pari
            if pari:
                l.append(pari)
            # print('pari', pari)
            l.extend(flatten_fingerprint(fp.copy(), pari.copy()))
            pari = []
    return l


if __name__ == "__main__":

    project = sys.argv[1]
    sample = sys.argv[2]

    file = os.path.join(project, f'sgs-prelim/mafft/{sample}.fasta')

    with open(file, 'r') as f:
        seqs = {record[0]: [char for char in record[1]] for record in SimpleFastaParser(f)}
    counts = get_counts(seqs)
    cons = determine_consensus(counts)
    # seqs = clean_pb_errors(seqs, counts, cons)

    # print(seqs)
    # print(len(seqs))
    for seq in seqs:
        umi = seq.split('_')[-2]

        fp_bin = os.path.join(project, f'bins/{sample}/aligned/{umi}.fasta')

        with open(fp_bin, 'r') as f:
            umi_seqs = {record[0]: [char for char in record[1]] for record in SimpleFastaParser(f)}

        print(umi, len(umi_seqs))
        fingerprints = get_pcr_fp(umi_seqs)
        print(fingerprints)

        # par = []
        # ffp = flatten_fingerprint(fingerprints, par)
        #
        # uniq_fp = []
        # for fp in ffp:
        #     if fp not in uniq_fp:
        #         uniq_fp.append(fp)
        #
        # for fp in uniq_fp:
        #     print(fp)


    # with open('tmp.fa', 'w') as f:
    #     for seq in seqs:
    #         seq_flat = ''.join(seqs[seq])
    #         f.write(f'>{seq}\n{seq_flat}\n')