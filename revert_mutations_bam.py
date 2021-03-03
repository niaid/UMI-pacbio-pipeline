#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 18:31:16 2020

@author: bayatmokhtarie

Note you need to convert your sam file to an indexed bam file.
Do so by

samtools view -S -b samfile_input.sam > bam_file_output.bam
samtools index bam_file_output.bam
samtools fasta test_revert_bam.bam > test_revert_bam.fasta
"""
#pylint: disable=C0413
import sys
import pickle
import itertools
import numpy as np
import pandas as pd
import pysam
from pysam import VariantFile

def load_no_change_dict(pickle_path, sample):
    """
    Loads the pickle dictionary and selects the sample from it

    Parameters
    ----------

    pickle_path: str
        Location of pickle dict

    sample: str
        The sample to load

    Returns
    -------

    no_change_dict: dict
        The dictionary to skip changes
        keys are positions, and variant, values are type
    """
    with open(pickle_path, 'rb') as handle:
        variant_dict = pickle.load(handle)
    if sample in variant_dict:
        no_change_dict = variant_dict[sample]
    else:
        no_change_dict = {}
    return no_change_dict

def make_no_change_dict(vcf_path):
    """
    Makes a no change list from a vcf file

    Parameters
    ----------
    vcf_path : str
        location to vcf file

    Returns
    -------
    no_change_list: list
        Locations to not change
    """
    vcf_in = VariantFile(vcf_path)
    no_change_dict = {}
    for rec in vcf_in.fetch():
        start = rec.pos
        ref = rec.ref
        key = tuple(range(start, start+len(rec.ref)))
        if rec.alts:
            for variant in rec.alts:
                if len(variant) > len(rec.ref):
                    var_type = "insert"
                else:
                    var_type = 'other'
                variant = variant + "-"*(len(ref) - len(variant))
                long_key = key + tuple([variant])
                no_change_dict[long_key] = var_type
        else:
            variant = ""
            variant = variant + "-"*(len(ref) - len(variant))
            var_type = 'other'
            long_key = key + tuple([variant])
            no_change_dict[long_key] = var_type
    return no_change_dict

def parse_cigar_seq(cigar_tuple, raw_read):
    """
    This generates an aligned sequnce with a cigar tuple and raw read

    Parameters
    ----------
    cigar_tuple : tuple of tuples
        cigar tuple
    raw_read : str
        raw_seqnece read

    Returns
    -------
    aligned_read: str
        aligned read
    """
    aligned_read = ""
    for key, count in cigar_tuple:
        if key == 0:
            aligned_read += raw_read[:count]
            raw_read = raw_read[count:]
        if key == 1:
            raw_read = raw_read[count:]
        if key == 2:
            aligned_read += '-'*count
        if key == 4:
            aligned_read += raw_read[:count]
            raw_read = raw_read[count:]
    return aligned_read

def parse_cigar_seq_inserts(cigar_tuple, raw_read):
    """
    This generates an aligned sequnce with a cigar tuple and raw read.
    It includes inserts

    Parameters
    ----------
    cigar_tuple : tuple of tuples
        cigar tuple
    raw_read : str
        raw_seqnece read

    Returns
    -------
    aligned_read: str
        aligned read
    """
    aligned_read = ""
    for key, count in cigar_tuple:
        if key == 0:
            aligned_read += raw_read[:count]
            raw_read = raw_read[count:]
        if key == 1:
            aligned_read += raw_read[:count]
            raw_read = raw_read[count:]
        if key == 2:
            aligned_read += '-'*count
        if key == 4:
            aligned_read += raw_read[:count]
            raw_read = raw_read[count:]
    return aligned_read

def make_raw_read(cigar_tuple, inserted_read, no_inserted_read):
    """
    Creates a raw read from the three pieces of info

    Parameters
    ----------
    cigar_tuple : list of tuples
        cigar sequence
    inserted_read : str
        aligned read including inserts
    no_inserted_read : str
        aligned read without inserts

    Returns
    -------
    raw_read : str
        raw read
    """
    raw_read = ""
    for key, count in cigar_tuple:
        if key in set({0, 2, 4}):
            raw_read += no_inserted_read[:count]
            inserted_read = inserted_read[count:]
            no_inserted_read = no_inserted_read[count:]
        if key == 1:
            raw_read += inserted_read[:count]
            inserted_read = inserted_read[count:]
    return raw_read

def find_inserts(old_cigar, insert_variants_list):
    """
    Finds locations of inserts to keep and inserts to remove

    Parameters
    ----------
    old_cigar : list of tuples
        cigar string
    no_change_dict : dict
        values given by vcf_file
    sam_file_ref_start : int
        offset of reference sequence

    Returns
    -------
    good_inserts : list
        list of inserts to keep
    bad_inserts : list
        List of inserts to remove
    """
    spot = 0
    inserts = []
    for key, count in old_cigar:
        if key == 1:
            inserts += list(range(spot, spot + count))
            spot += count
        else:
            spot += count
    good_inserts = [x for x in inserts if x in insert_variants_list]
    bad_inserts = [ x for x in inserts if x not in good_inserts]
    return good_inserts, bad_inserts

def build_new_cigar(read, good_inserts, bad_inserts):
    """
    Builds a new cigar sequence

    Parameters
    ----------
    read : str
        read with inserts and deletions
    insert_list : list
        location of inserts

    Returns
    -------
    cigar: list of tuples
        cigar sequence

    finished_read: str
        new read with insertions removed
    """
    finished_read = list(read)
    read = read.replace('-', '2')
    prep_read = list(read)
    for item in bad_inserts:
        prep_read[item] = '*'
        finished_read[item] = '*'
    for item in good_inserts:
        prep_read[item] = 1
    prep_read = [x for x in prep_read if x != '*']
    for num, item in enumerate(prep_read):
        if item != '2' and item != 1:
            prep_read[num] = 0
        if item == '2':
            prep_read[num] = int(item)
    finished_read = ''.join([x for x in finished_read if x != '*'])
    cigar = [(k, sum(1 for i in g)) for k, g in itertools.groupby(prep_read)]
    return cigar, finished_read

def revert_mutations(file_path, no_change_dict):
    """
    Reverts the mutations to the concensus in each position in change_list

    Parameters
    ----------
    samfile_path : str
        location of fasta file
    no_change_list : list
        list object with positions to revert mutations
    out_path: str
        place to write new changed fasta file
    Returns
    -------
    None.

    Side Effects
    ------------
    It writes a fasta file to the output location
    """
    #pylint: disable=E1101, R0914
    samfile = pysam.AlignmentFile(file_path, 'rb')
    sam_file_iter = samfile.fetch()
    sam_file_dict = {}
    new_sam_file_list = []
    sam_file_ref_start = []
    cigar_dict = {}
    pos_dict = {}
    for num, item in enumerate(sam_file_iter):
        cigar_tuple = item.cigartuples
        cigar_dict[item.query_name] = cigar_tuple
        sam_file_dict[item.query_name] = list(parse_cigar_seq(cigar_tuple,
                                                         item.query_sequence))
        new_sam_file_list.append(item)
        sam_file_ref_start += [item.reference_start]
        pos_dict[item.query_name] = item.pos

    # Change consencus columns first
    samfile_df = pd.DataFrame.from_dict(sam_file_dict, orient='index')
    mode_list = samfile_df.mode()
    col_num = len(list(samfile_df))

    sam_file_ref_start = sam_file_ref_start[0]
    no_change_list = []
    insert_variants_list = []
    for item in no_change_dict:
        if no_change_dict[item] == 'other':
            no_change_list += list(item[:-1])
        else:
            insert_length = len(item[-1])
            insert_variants_list += list(range(item[0], item[0] + insert_length))
    no_change_list = set(no_change_list)
    no_change_list = [x - sam_file_ref_start - 1 for x in no_change_list]
    insert_variants_list = set(insert_variants_list)
    insert_variants_list = [int(x) - int(sam_file_ref_start) for x in insert_variants_list]

    consensus = [x for x in range(col_num) if x not in no_change_list]
    con_mode = mode_list[consensus]
    samfile_df[consensus] = np.nan

    samfile_df[consensus] = samfile_df[consensus].fillna(con_mode.iloc[0])

    # next we need to find rows for variants,
    # then replace all other rows with consensus

    # REMBER TO SUBTRACT 1 because python starts at 0
    pairs = []
    for long_key, var_type in no_change_dict.items():
        if var_type == 'other':
            key = long_key[:-1]
            value = long_key[-1]
            key_list = [x - sam_file_ref_start - 1 for x in key]
            value = list(value)
            rows = samfile_df[key_list][(samfile_df[key_list] == value).all(1)].index
            pairs += list(itertools.product(rows, key_list))
    change_pairs = [x for x in itertools.product(samfile_df.index.tolist(), list(no_change_list))\
                    if x not in pairs]

    index_pos = dict(zip(samfile_df.index, range(len(samfile_df.index))))
    change_rows = [index_pos[x[0]] for x in change_pairs]
    change_cols = [x[1] for x in change_pairs]
    samfile_df.values[change_rows, change_cols] = np.nan
    samfile_df[no_change_list] = samfile_df[no_change_list].fillna(mode_list[no_change_list].iloc[0])
    samfile_df['read'] = samfile_df.sum(axis=1)

    new_reads = dict(zip(samfile_df.index,samfile_df['read']))
    new_cigars = {}
    for new_read in new_sam_file_list:
        header_id = new_read.query_name
        old_cigar = new_read.cigartuples
        if old_cigar[0][0] == 4:
            new_pos = new_read.pos - old_cigar[0][1]
            new_read.pos = new_pos
        good_inserts, bad_inserts = find_inserts(old_cigar, insert_variants_list)
        aligned_inserts =  parse_cigar_seq_inserts(old_cigar, new_read.query_sequence)
        aligned_no_inserts = new_reads[header_id]
        raw_read = make_raw_read(old_cigar, aligned_inserts, aligned_no_inserts)
        new_cigar, finished_read = build_new_cigar(raw_read, good_inserts, bad_inserts)
        new_cigars[header_id] = new_cigar
        new_read.cigartuples = new_cigar
        finished_read = finished_read.replace('-', '')
        new_read.query_sequence = finished_read
        if len(bad_inserts) > 0:
            print(f"{header_id}:{','.join([str(x) for x in bad_inserts])}", file=sys.stderr)

    outfile = pysam.AlignmentFile("-", "wb", template=samfile)
    for item in new_sam_file_list:
        outfile.write(item)

if __name__ == '__main__':
    FILE_PATH = sys.argv[1]
    DICT_PATH = sys.argv[2]
    SAMPLE = FILE_PATH.split('/')[-1].split('_')[0]
    if len(sys.argv) > 3 and sys.argv[3] == 'True':
        NO_CHANGE_DICT = make_no_change_dict(DICT_PATH)
    else: 
        NO_CHANGE_DICT = load_no_change_dict(DICT_PATH, SAMPLE)
    revert_mutations(FILE_PATH, NO_CHANGE_DICT)
