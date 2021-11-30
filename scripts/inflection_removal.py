#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function removes low count umis based on counts below the inflection point
"""
import sys
from itertools import groupby
import numpy as np
from Bio import SeqIO
from kneebow.rotor import Rotor
from umi_dedup import strip_umi

def find_kneeinf_cutoff(umi_counts):
    """
    returns the cutoff value corresponding to the inflection point

    Parameters:

    umi_counts: list
        list of umi counts sorted high to low

    Returns:

    cut_off: int that represents the cutoff
    """

    # rle run length encoding. list of tuples
    # First entry is umi count
    # second entry is frequency.
    if len(umi_counts) == 1:
        infl = umi_counts[0]
        knee = umi_counts[0]
        knee_fail = True
        return infl, knee, knee_fail
    rle = [(k, sum(1 for i in g)) for k, g in groupby(umi_counts)]
    rle.sort(key=lambda x: x[0], reverse=True)
    freqs = np.asarray([x[1] for x in rle])
    # dedup counts and sort
    umi_totals = list(set(umi_counts))
    umi_totals = np.asarray(sorted(umi_totals, reverse=True))
    run_rank = np.cumsum(freqs) - (freqs -1)/2
    log_rank = np.log10(run_rank)
    log_counts = np.log10(umi_totals)
    d1n = np.diff(log_counts)/np.diff(log_rank)
    # We used the estimated derivative to calculate the inflection point
    # Then count the count corresponding to the location of the inflection point
    right_edge = np.argmin(d1n)
    # Takes inflection point
    infl = umi_totals[right_edge]

    left_edge = np.argmax(d1n)
    log_rank = log_rank[left_edge:right_edge]
    log_counts = log_counts[left_edge:right_edge]
    data = np.column_stack((log_rank, log_counts))
    rotor = Rotor()
    if len(log_counts) > 0:
        rotor.fit_rotate(data)
        knee_idx = rotor.get_knee_index()
        knee = umi_totals[knee_idx + left_edge]
        knee_fail = False
    else:
        knee = infl
        knee_fail = True
    return infl, knee, knee_fail

def prep_headers(fasta_file, infl_yes):
    """
    This function takes a fasta file path and extracts the headers.
    It ranks the umis by count

    Parameters
    ----------
    fasta_path : list
        list of sequences in fasta file. each sequence is in Bio format

    infl_yes: boolean
        Boolean variable to us inflection point instead of knee
    Returns
    -------
    header_dict: dict
        wether or not to keep each header
    """

    headers = [[record.id] +  strip_umi(record.id).split('_') for record in fasta_file]
    headers.sort(key=lambda x: int(x[2]), reverse=True)
    header_dict = {}
    umi_counts = [int(x[2]) for x in headers]

    infl, knee, knee_fail = find_kneeinf_cutoff(umi_counts)
    for item in headers:
        if infl_yes:
            if int(item[2]) < infl:
                header_dict[item[0]] = False
            else:
                header_dict[item[0]] = True
        else:
            if int(item[2]) < knee:
                header_dict[item[0]] = False
            else:
                header_dict[item[0]] = True

    return header_dict, infl, knee, knee_fail

def main(project, file_path, umi_out_folder, infl_yes):
    """
    See module docstring

    Parameters
    ----------
    fasta_path : str
        path to fasta file

    Returns
    -------
    None.

    """
    if umi_out_folder == 'fake-umi-curation':
        fasta_path = project + umi_out_folder + "/" + file_path + "/" + file_path + "_post_nw.fasta"
    else:
        fasta_path = project + umi_out_folder + "/" + file_path + "/" + file_path + "_keep_nw.fasta"
    fasta_file = list(SeqIO.parse(fasta_path, "fasta"))
    print("Computing Inflection Point")
    header_dict, infl, knee, knee_fail = prep_headers(fasta_file, infl_yes)
    keeps = []
    remove = []
    for record in fasta_file:
        if header_dict[record.id]:
            keeps.append(record)
        else:
            remove.append(record)
    # Still need to write the files.
    # Write keeps and remove to appropraitely names files
    with open(project + umi_out_folder + "/" + file_path + "/" + file_path + "_final.fasta", "w") as handle:
        SeqIO.write(keeps, handle, 'fasta')
    with open(project + umi_out_folder + "/" + file_path + "/" + file_path + "_removes_inflection.fasta", "w") as handle:
        SeqIO.write(remove, handle, 'fasta')
    if knee_fail:
        knee_fail = "####knee point calculation failed####"
    else:
        knee_fail = ''
    with open(project + umi_out_folder + '/' + file_path + '/' + file_path + "infl_knee.txt", "w") as handle:
        print(f"{file_path},{knee},{infl}{knee_fail}", file=handle)

if __name__ == '__main__':
    PROJECT = sys.argv[1]
    FILE_PATH = sys.argv[2]
    UMI_OUT_FOLDER = sys.argv[3]
    if len(sys.argv) > 4 and sys.argv[4] == 'True':
        INFL = True
    else:
        INFL = False
    main(PROJECT, FILE_PATH, UMI_OUT_FOLDER, INFL)
