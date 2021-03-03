#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: bayatmokhtarie
"""

import sys
import subprocess
import os
from Bio import SeqIO


def main(project):
    """
    This module takes the fasta files and calculate read length
    Parameters
    ----------
    project : string
        path to the final read fasta files, there should be / at the end.
    file_path : string
        sample name without the extension.
    Returns
    -------
    None.
    """
    files = os.listdir(project + "finalreads")
    if not os.path.exists(project+"out_readlength"):
        subprocess.check_call(['mkdir', project+"out_readlength"])
    for file_path in files:
        path = project + "finalreads" + os.sep + file_path
        out_file = project + "out_readlength/"+ file_path.replace('read.fasta','_readlen.txt')
        lengths = map(len, SeqIO.parse(path, 'fasta'))
        with open(out_file, 'w') as read_count:
            for item in lengths:
                read_count.write(str(item) + '\n')

if __name__ == '__main__':
    PROJECT = sys.argv[1]
    main(PROJECT)
