#!/usr/bin/env python3
import functools
import json
import multiprocessing
import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import tqdm
import pysam
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import combinations
from seqlib import gen_coord_map, get_counts, determine_consensus


def parse_uc_file(fp_uc):
    clusters = dict()
    cluster_centroids = dict()
    with open(fp_uc, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] == 'C':
                cluster_id = int(line[1])
                read_id = line[8]
                cluster_centroids[cluster_id] = read_id
            else:
                if line[0] in ('S', 'H'):
                    cluster_id = int(line[1])
                    read_id = line[8]
                    if cluster_id in clusters:
                        clusters[cluster_id].append(read_id)
                    else:
                        clusters[cluster_id] = [read_id]
    return clusters, cluster_centroids


def assess_cluster_sizes(cluster_sizes):

    if not cluster_sizes:
        return None

    tot = sum([cs[1] for cs in cluster_sizes])
    cluster_sizes = sorted([(cs[0], cs[1]/tot) for cs in cluster_sizes], key=lambda cs: cs[1], reverse=True)

    if len(cluster_sizes) == 1:
        return cluster_sizes[0][0]
    else:
        if cluster_sizes[0][1] > 2*cluster_sizes[1][1] and cluster_sizes[0][1] >= 0.5 and cluster_sizes[0][1]*tot >= 10:
            return cluster_sizes[0][0]
    return None


if __name__ == "__main__":
    fp_uc = sys.argv[1]
    fp_ccs = sys.argv[2]
    fp_clusters_cons = sys.argv[3]
    fp_cc_cons = sys.argv[4]
    fp_out = sys.argv[5]

    umi = fp_uc.split('/')[-1].strip('.uc')
    clusters, centroids = parse_uc_file(fp_uc)
    core_cluster = assess_cluster_sizes([(c, len(clusters[c])) for c in clusters])

    if core_cluster is not None:

        with open(fp_ccs, 'r') as f:
            ccs_seqs = {record[0]: record[1] for record in SimpleFastaParser(f)}

        with open(fp_out, 'w') as f:
            for read in clusters[core_cluster]:
                f.write(f'>{read}\n{ccs_seqs[read]}\n')

        with open(fp_clusters_cons, 'r') as f:
            cons_seqs = {record[0]: record[1] for record in SimpleFastaParser(f)}

        cons_read_id = None
        for read in cons_seqs:
            if read.__contains__(centroids[core_cluster]):
                cons_read_id = read
        if not cons_read_id:
            print(f'Error obtaining consensus sequence for: {umi}')
            sys.exit()

        cc_seq = cons_seqs[cons_read_id]

        with open(fp_cc_cons, 'w') as f:
            f.write(f'>{umi}-CC{core_cluster}-cons\n{cc_seq}\n')



