"""
This script uses header file to subset the sequences and write them as fasta file

project: top level directory of the project, must include "/"
file_path: fasta file to be subsampled.
wanted_file: txt file contains headers we want to take sequences from
result_file: fasta file that is subsampled from file_path using headers from
wanted_file.

"""

import sys
from Bio import SeqIO

def main(file_path, wanted_file,result_file):
    """
    This module takes the  fasta and txt files and return
    fasta file.

    Parameters
    ----------
    file_path and wanted_file : string
        path to the fasta and txt files. There should be / at the end.
    file_name : string
        concat. fasta file without extension.
    Returns
    -------
    None.

    """
    fasta_file = list(SeqIO.parse(file_path, "fasta"))
    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)
    keeps = []
    for seq in fasta_file:
        if any(item in seq.id for item in wanted):
            keeps.append(seq)
    with open(result_file, "w") as f:
        SeqIO.write(keeps, f, "fasta")



if __name__ == '__main__':
    FILE_PATH = sys.argv[1]
    WANTED_FILE = sys.argv[2]
    RESULT_FILE = sys.argv[3]
    main(FILE_PATH,WANTED_FILE, RESULT_FILE)

