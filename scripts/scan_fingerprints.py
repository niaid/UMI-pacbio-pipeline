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


def clean_pb_errors(seq_dict, counts, cons):
    for index, counter in enumerate(counts):
        for char in counter:
            if counter[char] == 1 or counter[char] < 0.01*sum(counter.values()):
                for seq in seq_dict:
                    if seq_dict[seq][index] == char and index in cons:
                        seq_dict[seq][index] = cons[index]
    return seq_dict


def get_largest_vars(seq_dict):
    counts = get_counts(seq_dict)
    putative_vars = []
    for index, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c < 0.9 and counter.most_common(2)[1][0] not in '-' and counter.most_common(2)[0][0] not in '-':
            if counter.most_common(2)[1][1] >= 3:
                putative_vars.append((index, counter.most_common(2)[1][0],
                                      counter.most_common(2)[1][1] / sum(counter.values()),
                                      counter.most_common(2)[0][0]))

    if putative_vars:
        max_counts = max(c[2] for c in putative_vars)
        # print(max_counts)
    else:
        return []
    vs = list(filter(lambda c: c[2] > 0.9*max_counts, putative_vars))
    return vs


def get_seqs_w_fp(seq_dict, ls):
    s = dict()
    sn = dict()
    for seq in seq_dict:
        match = True
        for i, site in enumerate(ls):
            #print(seq, site[0], seq_dict[seq][site[0]], site[1])
            if seq_dict[seq][site[0]] != site[1]:
                match = False
                break
        if match:
            s[seq] = seq_dict[seq]
        else:
            sn[seq] = seq_dict[seq]
    return s, sn


def assess_linkage(sites, seq_dict):
    if not sites:
        return []
    linkage_map = {seq: [0] * len(sites) for seq in seq_dict}
    for seq in seq_dict:
        for i, site in enumerate(sites):
            if seq_dict[seq][site[0]] == site[1]:
                linkage_map[seq][i] = 1
    linkage_map = np.array([linkage_map[seq] for seq in linkage_map])
    lsites_indices = assess_linkage_map(linkage_map)
    # if (1008 in (site[0] for site in sites)):
    #     print(lsites_indices)
    #     plt.imshow(linkage_map)
    #     plt.show()
    return [[sites[i] for i in cc] for cc in lsites_indices]


def compute_link(a, b):
    return sum(ai*bi == 1 for ai, bi in zip(a, b))/sum(ai or bi for ai, bi in zip(a, b))


def assess_linkage_map(lm):
    ns = lm.shape[1]
    nodes = list(range(ns))

    G = nx.Graph()
    G.add_nodes_from(nodes)

    for combo in combinations(list(range(ns)), 2):
        if compute_link(lm[:, combo[0]], lm[:, combo[1]]) > 0.9:
            G.add_edge(combo[0], combo[1])

    linked_indices = list(nx.connected_components(G))
    return linked_indices


def get_pcr_fp(seq_dict, round=0):

    if len(seq_dict) <= 4:
        return []
    fp = []
    sites = get_largest_vars(seq_dict)
    if not sites:
        return []

    lsites = assess_linkage(sites, seq_dict)[0]
    if len(lsites) > 4:
        return []
    # print(round, len(seq_dict), lsites)
    fp.extend([lsites])

    lseqs, lseqs_n = get_seqs_w_fp(seq_dict.copy(), lsites)
    # if len(lseqs) == 0:
    #     return []
    fp.extend([lsites + sites for sites in get_pcr_fp(lseqs, round=round+1)])

    # print(round, fp)
    if round == 0:
        if lsites and lsites[0][2] > 0.3:
            rlsites = [(ls[0], ls[3], 1-ls[2], ls[1]) for ls in lsites]
            fp.extend([rlsites] + get_pcr_fp(lseqs_n, round=round+1))
    else:
        fp.extend(get_pcr_fp(lseqs_n, round=round+1))
    return fp


def load_sam_seqs(fp_sam):

    samfile = pysam.AlignmentFile(fp_sam, "r")

    samfile.check_index()
    ri = samfile.get_reference_length(samfile.get_reference_name(0))

    sam_seqs = dict()
    try:
        for read in samfile.fetch():

            ms = ''
            ap = read.get_aligned_pairs()

            aln_dict = {a: None for a in range(ri)}

            for pp in ap:
                if pp[1] or pp[1] == 0:
                    aln_dict[pp[1]] = pp[0] if pp[0] else None

            for i in aln_dict:
                qp = aln_dict[i]
                if qp:
                    ms += read.seq[qp]
                else:
                    ms += '-'

            sam_seqs[read.qname] = [a for a in ms]
    except StopIteration as e:
        print(f'Warning: problem with umi {fp_sam}')
        print(e)
    return sam_seqs


def load_cons_seq(fp):
    with open(fp, 'r') as f:
        cseq = [c.upper() for c in next(SimpleFastaParser(f))[1]]
        # print(cseq)
    return cseq


def fingerprint(umi_data, project, sample, sgs_cons):
    umi_label = umi_data[0]
    umi = umi_label.split('_')[-2]
    umi_sgs_aligned = umi_data[1]
    umi_sgs_aligned = [a.upper() for a in umi_sgs_aligned]
    # fp_bin = os.path.join(project, f'bins/{sample}/aligned/{umi}.fasta')
    # with open(fp_bin, 'r') as f:
    #     umi_seqs = {record[0]: [char for char in record[1]] for record in SimpleFastaParser(f)}

    umi_seqs = load_sam_seqs(os.path.join(project, f'bins/{sample}/aligned/{umi}-sorted.sam'))

    if not umi_seqs:
        fingerprints = []
    else:
        umi_counts = get_counts(umi_seqs)
        umi_cons = determine_consensus(umi_counts, threshold=0.9)
        umi_seqs = clean_pb_errors(umi_seqs, umi_counts, umi_cons)
        try:
            fingerprints = get_pcr_fp(umi_seqs)
        except Exception:
            print(umi, 'Recursion error.')
            raise

    dd = []
    for i, s in enumerate(umi_sgs_aligned):
        cb = i in sgs_cons
        if cb and s != sgs_cons[i]:
            dd.append((i, s))

    # local_seq = ''.join(determine_consensus(umi_counts, threshold=0).values())
    local_seq = load_cons_seq(os.path.join(project, f'bins/{sample}/consensus/{umi}.fasta'))

    glmap_dict = gen_coord_map(umi_sgs_aligned, local_seq)
    lgmap_dict = gen_coord_map(local_seq, umi_sgs_aligned)

    return {umi: {'fp_pcr': fingerprints, 'fp_dd': dd}}, {umi: glmap_dict}, {umi: lgmap_dict}


if __name__ == "__main__":

    project = sys.argv[1]
    sample = sys.argv[2]
    if len(sys.argv) == 4:
        threads = int(sys.argv[3])
    else:
        threads = multiprocessing.cpu_count()

    sgs_file = os.path.join(project, f'sgs-prelim/mafft/{sample}.fasta')

    with open(sgs_file, 'r') as f:
        sgs_seqs = {'_'.join(record[0].split('_')[-2:]): [char for char in record[1]] for record in SimpleFastaParser(f)}
    sgs_counts = get_counts(sgs_seqs)
    sgs_cons = determine_consensus(sgs_counts)
    sgs_cons = {i: sgs_cons[i].upper() for i in sgs_cons}

    with open(os.path.join(project, f'sgs-prelim/sgs-consensus/{sample}.fasta'), 'w') as f:
        f.write(f'>{sample}_sgs-consensus\n')
        for i in range(len(next(iter(sgs_seqs.values())))):
            if i in sgs_cons:
                f.write(sgs_cons[i])
            else:
                f.write(r'X')

    umi_data = [(umi_label, sgs_seqs[umi_label]) for umi_label in sgs_seqs.keys()]
    # print([(i, umi[0]) for i, umi in enumerate(umi_data)])

    pool = multiprocessing.Pool(processes=threads)

    pcr_fp_dict = {}
    gl_map_dict = {}
    lg_map_dict = {}

    with tqdm.tqdm(total=len(umi_data), unit='UMI') as pbar:
        worker = pool.imap_unordered(functools.partial(fingerprint, project=project, sample=sample, sgs_cons=sgs_cons),
                                     sorted(umi_data, key=lambda ud: int(ud[0].split('_')[1]), reverse=True))
        try:
            for result in worker:
                # print(result[0].keys())
                pcr_fp_dict.update(result[0])  # Update fingerprints
                gl_map_dict.update(result[1])  # Update fingerprints
                lg_map_dict.update(result[2])  # Update fingerprints
                pbar.update()
            pool.close()
            pool.join()
        except Exception as e:
            pool.terminate()
            raise

    with open(os.path.join(project, f'sgs-prelim/fingerprints/{sample}.json'), 'w') as f:
        json.dump([pcr_fp_dict, gl_map_dict, lg_map_dict], f, indent=2)
