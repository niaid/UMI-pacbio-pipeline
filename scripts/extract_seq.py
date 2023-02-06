#!/usr/bin/env python
"""
This module collects all the sequences with the same UMI
then places them in a file together.
"""
# -*- coding: utf-8 -*-


import sys
import os
import subprocess
from Bio import SeqIO

#Extracts reads based on associated UMI sequence#
def main(project, file_path, number_file):
    """ See module doc string """
    fasta_file = project + "umi_stats/" + file_path + "_UMI.fasta"  # Input fasta file
    result_file = "seq.fasta" # Output fasta file

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
    end = False

    with open(project +"sequences/" + file_path + "/" + \
              str(number_file)+"_"+result_file, "w") as eq_umi_file:
        for seq in fasta_sequences:
            seq_id = seq.id.split("_")
            umi = seq_id[3]
            if umi == number_file:
                SeqIO.write([seq], eq_umi_file, "fasta")


if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    NUMBER_FILE = sys.argv[3]
    print(NUMBER_FILE)
    main(PROJECT, FILE_PATH, NUMBER_FILE)
