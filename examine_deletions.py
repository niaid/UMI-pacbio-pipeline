#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates a table to detect inframe deletions,
multiple 1-base deletions,
"""

import sys
import pysam
from revert_mutations_bam import parse_cigar_seq_inserts

def main(bam_file_input, outpath):
    """
    This function reads the bam file input and tells us about inframe/out of
    frame deletions

    Parameters
    ----------
    bam_file_input : str
        Location of bam file

    Returns
    -------
    None.

    """
    bamfile = pysam.AlignmentFile(bam_file_input, 'rb')
    bam_file_iter = bamfile.fetch()
    list_to_write = []
    hap_dict = {}
    for record in bam_file_iter:
        cigar_tuple = record.cigartuples
        deletions = [x for x in cigar_tuple if x[0] == 2]
        bad_deletions = deletions.count((2, 1))
        seq_read = parse_cigar_seq_inserts(cigar_tuple, record.query_sequence)
        if seq_read in hap_dict:
            hap_dict[seq_read] += 1
        else:
            hap_dict[seq_read] = 1
        for _, count in deletions:
            if count > 2:
                inframe = count %3 == 0
                list_to_write.append((record.query_name,\
                                      inframe, count, bad_deletions, seq_read))
    uni_hap_dict = {}
    for num, item in enumerate(hap_dict.items()):
        hap = item[0]
        uni_hap_dict[hap] = "H" + str(num)
    with open(outpath, "w") as handle:
        print(f"header,inframe,deletion length,#of 1-base deletions,hap group", file=handle)
        for item in list_to_write:
            item = [str(x) for x in item]
            print(f"{','.join(item[:-1])},{uni_hap_dict[item[-1]]}", file=handle)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
