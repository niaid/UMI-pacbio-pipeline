#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module find haplotypes (unique sequences) and find the frequency of each
of them. The results are written as fasta format where the headers are sample
name+ frequency of occurance.
"""
import sys
import os
import itertools

from Bio import SeqIO


def main(project, write_folder):
    """
    This module takes the fasta files and return fasta files where
    in headrs are the sample name + frequency of each sequence and the
    sequence as read.
    Parameters
    ----------
    project : string
        Path where the fasta files are.
    write_folder : string
        path to the write files.
    Returns
    -------
    None.

    """

    files = os.listdir(project)
    files=[x for x in files if '.fasta' in x]

    for sample in files:
        fasta_file = list(SeqIO.parse(project + sample, "fasta"))
        sam_id = sample.split("_")[0]
        hap_dict = {}
        hap_counts = {}
        for record in fasta_file:
            if record.seq in hap_counts:
                hap_counts[record.seq] += 1
            else:
                hap_counts[record.seq] = 1
        ones = [x for x in hap_counts.items() if x[1] == 1]
        not_one = [x for x in hap_counts.items() if x[1] != 1]
        not_one.sort(key=lambda x: x[1])
        total_haps = not_one + ones
        hap_dict = dict(zip([x[0] for x in total_haps],\
                              [f"H{x}" for x in range(1,len(total_haps) + 1)]))
        with open(write_folder + '/' + sam_id + "hap-groups.csv", "w") as handle:
            print("header,hap-group", file=handle)
            for record in fasta_file:
                print(f"{record.id},{hap_dict[record.seq]}", file=handle)

        sequences = [str(x.seq) for x in fasta_file]
        sequences.sort()
        counts = [(k, sum(1 for i in g)) for k, g in itertools.groupby(sequences)]
        with open(write_folder + "/" + sam_id + '_uniqueseq.fasta', 'w') as outfile:
            for seq, count in counts:
                    print(f">{sam_id}-count-{hap_dict[seq]}:{count}"+"\n"+seq, file=outfile)

if __name__ == '__main__':
    PROJECT = sys.argv[1]
    WRITE_FOLDER = sys.argv[2]
    main(PROJECT, WRITE_FOLDER)

#   project = '/Users/bayatmokhtarie/Desktop/consensus_revert/AMOEBA/consensus_revert/reverted_mutations/mafft/'
