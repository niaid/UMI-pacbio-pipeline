#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This extracts reads with UMIS that appear more than 3 times.
"""

#pylint: disable=R1704, R0914, C0103, E0401
import sys
import subprocess
import csv
import os
from inflection_removal import find_kneeinf_cutoff


def getPrimerId(trimmed_FP_RP_file, primerId_file, RT_primer_RC):
    """
    Function reads the trimmed files and looks for 8 base primerID sequence after
    RT primer. Each primerID is written to a csv file with associated sequence ID
    and sequence

    Parameters
    ----------
    trimmed_FP_RP_file : str
        path to trimmed reads file. It has no Forward or reverse primers
    primerId_file : str
        path to csv file that stores UMI, header, read
    RT_primer_RC : str
        reverse complement of RT primer.

    Returns
    -------
    None.

    """
    myTup = ()
    with open(trimmed_FP_RP_file, 'r') as trimmed, \
        open(primerId_file, 'w') as UMI_csv_file:
        for line in trimmed:
            if line.startswith('>'):
                name = line.strip()
            else:
                value = line
                detected_RTPrimer_RC = value[len(value)-(len(RT_primer_RC) \
                                                 + 9):len(value)-9]
                if RT_primer_RC == detected_RTPrimer_RC:
                    UMI = value[len(value)- 9 :len(value)]
                    myTup = (UMI.strip(), value.strip())
                    writer = csv.writer(UMI_csv_file)
                    writer.writerow([myTup[0], name, myTup[1]])



def countPrimerId(primerId_file, count_pid):
    """
    Function reads the primerID file generated in getPrimerId() and associates a
    count column for each occurence of a primerID sequence
    """

    counter = {}
    data = []
    with open(primerId_file, 'r') as primerId_file:
        for row in csv.reader(primerId_file):
            key = row[0]
            if key not in counter:
                counter[key] = 1
            else:
                counter[key] += 1
            row.insert(0, counter[key])
            data.append(row)
    for row in sorted(data, key=lambda x: (x[1], x[0])):
        with open(count_pid, 'w') as count_pid_csv:
            fileWriter = csv.writer(count_pid_csv)
            fileWriter.writerow(row)



def csvToFasta(UMI_csv_file, pid_fasta):
    """
    Function reads the primerID file and writes to a fasta format file in the
    format seqID_UMI and then the sequence

    Parameters
    ----------
    UMI_csv_file : str
        path to the csv file with UMI, name, read
    pid_fasta : str
        fasta output of the csv file

    Returns
    -------
    None.

    """
    with open(UMI_csv_file, 'r') as pid_csv, open(pid_fasta, 'w') as pid_fasta:
        for line in pid_csv:
            value = line.split(",")
            pid_fasta.write(value[1]+'_'+value[0])
            pid_fasta.write('\n'+value[2])



def pid(pid_fasta, pidFile):
    """
    Function reads the fasta file generated from csvToFasta and outputs a file with
    list of UMI's

    Parameters
    ----------
    pid_fasta : str
        fasta output of the csv file
    pidFile : str
        file with all the unique UMIs

    Returns
    -------
    None.

    """
    with open(pid_fasta, 'r') as pid_fasta, open(pidFile, 'w') as pidFile:
        for line in pid_fasta:
            if line.startswith(">"):
                value = line.split("_")
                pid_value = value[3]
                pidFile.write(pid_value)


def umiSeq(pidcounts, umiSequences):
    """
    Read the counts umi file and write umi sequences to a file and 
    this step discards umi's <= inflection point

    Parameters
    ----------
    pidcounts : str
        path to UMI counts table.
    umiSequences : str
        Path to write UMIs that > inflection point.

    Returns
    -------
    None.

    """
    UMI_count_list = []
    with open(pidcounts,'r') as pid_counts_file:
        for line in pid_counts_file:
            umi_count = int(line.strip().split()[0])
            UMI_count_list.append(umi_count)
    print(UMI_count_list)
    if len(UMI_count_list) > 1 and max(UMI_count_list) > 3:
        cutoff, knee, knee_fail = find_kneeinf_cutoff(UMI_count_list)
        print(cutoff)
    
        with open(pidcounts, 'r') as pid_counts_file, open(umiSequences, 'w') as output:
            for line in pid_counts_file:
                line = line.strip(" ")
                if int(line.strip().split()[0]) >= cutoff:
                    seq = line.split(" ")
                    output.write(seq[1])
    else:
        cutoff = None
        knee = None
        knee_fail = None
    return cutoff, knee, knee_fail


def main(RT_primer_RC, project, file_path):
    """ See module doc string"""
    trimmed_FP_RP_file = project + "/trimmed/" + file_path \
      + "_trimmed.fasta"  # Input fasta file
    if not os.path.exists(project + "/umi_stats"):
        subprocess.check_call(['mkdir', project + "/umi_stats/"])
    if not os.path.exists(project + "/sequences"):
        subprocess.check_call(['mkdir', project + "/sequences/"])
    # Prepare output variables
    primerId_file = project + "/umi_stats/" + file_path + "_UMI_seq.csv"
    #count_pid = project + "/umi_stats/" + file_path + "_UMI_count.csv"
    pid_fasta = project + "/umi_stats/" + file_path + "_UMI.fasta"
    pidFile = project + "/umi_stats/" + file_path + "_UMI.txt"
    pidcounts = project + "/umi_stats/" + file_path + "_counts_UMI.txt"
    umiSequences = project + "/umi_stats/" + file_path + "_umi_seq.txt"
    #result_file = "seq.fasta" # Output fasta file

    getPrimerId(trimmed_FP_RP_file, primerId_file, RT_primer_RC)
    csvToFasta(primerId_file, pid_fasta)
    pid(pid_fasta, pidFile)

    #sort and unique the UMI to get occurence of UMI

    pidcounts_output = open(pidcounts, 'w+')
    sort_pid_arg = ['sort', pidFile]
    unique_pid_arg = ['uniq', '-c']
    sort_unique_pid_arg = ['sort', '-k', '1n']
    sort_pid = subprocess.Popen(sort_pid_arg, \
        stdout=subprocess.PIPE, shell=False)
    unique_pid = subprocess.Popen(unique_pid_arg, stdin=sort_pid.stdout, \
        stdout=subprocess.PIPE, shell=False)
    sort_unique_pid = subprocess.Popen(sort_unique_pid_arg, \
        stdin=unique_pid.stdout, stdout=pidcounts_output, shell=False)
    sort_pid.stdout.close()
    unique_pid.stdout.close()
    pidcounts_output.close()
    sort_unique_pid.communicate()[0]

    cutoff_file = project + "/umi_stats/inflection_points/" + file_path \
      + "_infl.txt"
    if not os.path.exists(project + "/umi_stats/inflection_points"):
        subprocess.check_call(['mkdir', project + "/umi_stats/inflection_points"])
    inflection, knee, knee_fail = umiSeq(pidcounts, umiSequences)
    if knee_fail:
        knee_fail = "###knee point failed###"
    else:
        knee_fail = ""
    with open(cutoff_file,'w') as f:
        f.write(file_path + "," + str(knee) + "," + str(inflection)+knee_fail)

if __name__ == '__main__':
    RT_PRIMER_RC = sys.argv[1]
    PROJECT = sys.argv[2]
    FILE_PATH = sys.argv[3]
    main(RT_PRIMER_RC, PROJECT, FILE_PATH)
