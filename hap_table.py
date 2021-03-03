#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This groups haplotypes across longitudinal samples
"""


import sys
import os


from Bio import SeqIO

def main(project, outpath):
    """
    Groups haplotypes for a sample

    Parameters
    ----------
    project : str
        location of samples
    outpath : str
        place to write file

    Returns
    -------
    None.
    """
    if os.path.isdir(project):
        files = os.listdir(project)
        files=[x for x in files if '.fasta' in x]
    else:
        files = [project.split('/')[-1]]
        project = '/'.join(project.split('/')[:-1]) + '/'
    hap_set = set()
    with open(project + outpath, 'w') as handle:
        print('Sample,Count,Perc,Haplotypes', file=handle)
    to_be_written = []
    for sample in files:
        fasta_file = list(SeqIO.parse(project + sample, "fasta"))
        sam_id = sample.split("_")[0]
        hap_counts = {}
        for record in fasta_file:
            if record.seq not in hap_set:
                hap_set.add(record.seq)
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
        total = sum(hap_counts.values())
        for seq, count in hap_counts.items():
            to_be_written.append([sam_id, f"{count}", f"{round(count/total,3)}", seq])
    hap_dict = dict(zip(list(hap_set), [ f"H{i}" for i in range(1,1+len(hap_set))]))
    with open(project + outpath, 'a') as handle:
        for item in to_be_written:
            print(f"{item[0]},{item[1]},{item[2]},{hap_dict[item[3]]}", file=handle)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])