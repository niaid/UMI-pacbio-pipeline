"""
This module selects the consensus read for each UMI with the highest read abundance
It write all the consensus reads into one fasta file.
"""

import sys
import os
import shutil
import subprocess
import numpy as np
from Bio import SeqIO

#pylint: disable=R0914
def read_cluster_size(cluster_stats_file, consensus_file, cluster_stats_out_file, umi):
    """
    This reads a .uc file and finds the consensus read with the most reads.
    It writes the consensus reads with the abundance counts.

    Parameters
    __________

    cluster_stats_file: str
        the path to the .uc file
    consensus_file: str
        path to the consensus reads file
    cluster_stats_out_file: str
        The path where the counts will be written
    umi: str
        the UMI

    Returns
    _______

    None
    """
    id_list = []
    unique_counts_list = []
    cluster_number_list = []
    with open(cluster_stats_file, 'r') as fin, open(consensus_file, 'r') as \
      fin2, open(cluster_stats_out_file, 'a') as fout:
        for line in fin:
            if line.startswith('C'):
                value = line.strip().split('\t')
                cluster_number = value[1]
                unique_counts = value[2]
                header = value[8]
                cluster_number_list.append(int(cluster_number))
                unique_counts_list.append(int(unique_counts))
                id_list.append(header)
        total = sum(unique_counts_list)
        current_line = 0
        for line in fin2:
            if line.startswith('>Cluster') and str(cluster_number_list[current_line]) in line:
                fout.write(line.strip() + "_" + str(total)+ "_" + umi + "_" + \
                           str(unique_counts_list[current_line]) +'\n')
                current_line += 1
            else:
                fout.write(line.strip() +'\n')

def main(project, file_path, number_file):
    """ See module doc string"""
    tau1 = .5 # Ratio of cluster reads to total UMI reads.
              # Value greater or equal .5 means at least 50% of the reads are in the first cluster
    tau2 = .5 # Ratio of second largest cluster to first largest cluster.
              # Value less than or equal to .5 means the first cluster was at least twice as big as the second
    tau3 = .95 # Value to assess possible UMI collision. Saying top clusters had a very similar number

    if not os.path.exists(project + "sequences/" + \
                          file_path + "/cluster_stats"):
        subprocess.check_call(['mkdir', project + "sequences/" + \
                               file_path + "/cluster_stats"])
        subprocess.check_call(['mkdir', project + "sequences/" + \
                               file_path + "/cluster_with_read_counts"])
    if not os.path.exists(project + "umi_collision"):
        subprocess.check_call(['mkdir', project + "umi_collision"])

    if not os.path.exists(project + 'error_insert'):
        subprocess.check_call(['mkdir', project + "error_insert"])
    
    cluster_file_path = project + "sequences/" + file_path
    fasta_file = cluster_file_path + "/cluster_stats/"+ str(number_file)+"_read.fasta"
    result_file = cluster_file_path + "/cluster_with_read_counts/" + \
    str(file_path) +"read.fasta"
    umi_collision = project + "umi_collision/" + str(file_path) + "collision.txt"
    error_insert =  project + "error_insert/" + str(file_path) + "error_insert.txt"
    umi_collision_reads = project + "umi_collision/" +str(file_path) + "collision.fasta"
    error_insert_reads = project + "error_insert/" + str(file_path) + "error_insert.fasta"
    id_list = []
    read_count_list = []
    seq_list = []
    list_index = 0

    read_cluster_size(cluster_file_path + "/cluster_usearch/" + number_file +\
                "_cluster.uc", cluster_file_path + "/consensus_usearch/" +\
                number_file + "_consensus.fasta", cluster_file_path + \
                "/cluster_stats/" + number_file + "_read.fasta", number_file)

    fasta_sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    with open(fasta_file, 'r') as fin, open(result_file, 'a') as fout,\
    open(umi_collision, 'a') as umi_collision_file, open(error_insert, 'a') as error_file,\
    open(umi_collision_reads, 'a') as umi_col_fasta, open(error_insert_reads, 'a') as error_fasta:
        for line in fin:
            if line.startswith('>'):
                id_list.append(line)
                header_id = line.split("_")
                read_count = header_id[3]
                read_count_list.append(int(read_count))
            else:
                seq_list.append(line)

        total = sum(read_count_list)
        list_index = np.argmax(read_count_list)
        str_read_count  = [str(x) for x in read_count_list]
        # We should add a bit of code here to select further using thresholds
        if float(read_count_list[list_index])/total >= tau1:  
            fout.write('>' + file_path + '_' + fasta_sequences[list_index].id)
            fout.write('\n'+str(fasta_sequences[list_index].seq)+'\n')
        else:
            sorted_reads = sorted(read_count_list, reverse=True)
            cluster_ratio = float(sorted_reads[1])/float(sorted_reads[0])
            if  cluster_ratio <= tau2:
                fout.write('>' + file_path + '_' + fasta_sequences[list_index].id)
                fout.write('\n'+str(fasta_sequences[list_index].seq)+'\n')
            elif cluster_ratio >= tau3:
                umi_collision_file.write(number_file + ',' + str(total) + "," + ",".join(str_read_count) + '\n')
                umi_col_fasta.write('>' + file_path + '_' + fasta_sequences[list_index].id)
                umi_col_fasta.write('\n'+str(fasta_sequences[list_index].seq)+'\n')
                shutil.copy(cluster_file_path + "/" + number_file + "_seq.fasta",\
                    project + "umi_collision/" + str(file_path) + "_" + number_file +".fasta")
            else:
                error_file.write(number_file + ',' + str(total) + "," + ",".join(str_read_count) + '\n')
                error_fasta.write('>' + file_path + '_' + fasta_sequences[list_index].id)
                error_fasta.write('\n'+str(fasta_sequences[list_index].seq)+'\n')
                shutil.copy(cluster_file_path + "/" + number_file + "_seq.fasta",\
                    project + "error_insert/" + str(file_path) + "_" + number_file + ".fasta")
if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    NUMBER_FILE = sys.argv[3]
    main(PROJECT, FILE_PATH, NUMBER_FILE)
