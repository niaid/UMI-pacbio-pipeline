"""
This script separates reads from possible umi collision into specific clusters

Inputs are

project: top level directory of the project, must include "/"
file_path: the sample id being examined. no extension
umi: this is the umi to copy from. Read from standard input using xargs
"""

import sys
import os
import subprocess
from Bio import SeqIO

def main(project, file_path, umi):
    """
    See module doc strings
    """
    uc_file = project + "sequences/" + file_path + "/cluster_usearch/" + umi + "_cluster.uc"

    header_dict = {}
    with open(uc_file, 'r') as handle:
        for line in handle:
            if line[0] != "C":
                line = line.split("\t")
                header_dict[line[8]] = line[1]

    # makes large folder for cluster reads
    cluster_reads_folder = project + "umi_collision/" + file_path + "/" + umi + "_cluster_reads/"
    if not os.path.exists(project + "umi_collision/" + file_path + "/"):
        subprocess.check_call(['mkdir', project + "umi_collision/" + file_path + "/"])
    if not os.path.exists(cluster_reads_folder):
        subprocess.check_call(['mkdir', cluster_reads_folder])
    cluster_dict = {}
    fasta_sequences = SeqIO.parse(project + "sequences/" + file_path + "/" + umi + "_seq.fasta", 'fasta')
    for record in fasta_sequences:
        cluster = header_dict[record.id]
        if cluster in cluster_dict:
            cluster_dict[cluster] += [record]
        else:
            cluster_dict[cluster] = [record]
    for cluster_num in cluster_dict:
        SeqIO.write(cluster_dict[cluster_num], cluster_reads_folder + "Cluster" + str(cluster_num) + ".fasta", "fasta")

if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    UMI = sys.argv[3]
    main(PROJECT, FILE_PATH, UMI)
