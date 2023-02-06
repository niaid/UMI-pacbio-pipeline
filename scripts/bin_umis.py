"""
bin_umis.py input rt_primer output_dir umi_length

Parses a file of oriented and PCR-primer-trimmed reads, where the reads thus end with a UMI sequence of length
"umi_length" preceded by the RT primer sequence "rt_primer".

UMIs are binned and all bins with size greater than 10 CCS are written to a respective FASTA file.
"""
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from scipy.ndimage import gaussian_filter1d
from seqlib import hamming

matplotlib.use('Agg')


def bin_umis(seqs, rtp_sequence, umi_len, fp_umi_labeled, qc=True):

    umi_map_dict = dict()
    read_count = 0
    exceptions = []

    with open(fp_umi_labeled, 'w') as f:
        for name, seq in seqs:
            read_count += 1
            rtp_detected = seq[-(len(rtp_sequence) + umi_len):-umi_len]  # Get RT primer sequence
            if hamming(rtp_detected, rtp_sequence) <= 2:  # Check if consistent with designed RTP sequence
                umi = seq[-umi_len:]
                if umi in umi_map_dict:
                    umi_map_dict[umi].append(name)
                else:
                    umi_map_dict[umi] = [name]

                f.write(f'>{name}_{umi}\n{seq}\n')
            else:
                exceptions.append(name)  # Track reads without valid RTP sequence

    if qc:
        print("\tNumber of reads with invalid RT primer sequence: {} ({:.2f}%)".format(
            len(exceptions),
            100*len(exceptions) / read_count))

    return umi_map_dict


def main(fp_in, rtp_rc, output_folder, umi_length):

    with open(fp_in, 'r') as fp_in_handler:
        seqs = SimpleFastaParser(fp_in_handler)
        filename = fp_in.split('/')[-1]
        sample_name = filename.replace('.fasta', '').replace('_trimmed', '')
        umi_map = bin_umis(seqs, rtp_rc, umi_length, os.path.join(output_folder, f'{sample_name}_UMI.fasta'))

    sorted_umis = [(len(umi_map[umi]), umi) for umi in sorted(umi_map, key=lambda u: (len(umi_map[u]), u))]

    with open(os.path.join(output_folder, f'{sample_name}_counts.txt'), 'w') as f_out:
        for count, umi in sorted_umis:
            f_out.write(f'{count}\t{umi}\n')

    counts = np.array([count for count, _ in reversed(sorted_umis)], dtype=float)
    smoothed_counts = gaussian_filter1d(counts, sigma=len(counts)/25)

    filter_threshold = 10

    plt.plot(counts, lw=2, color='b', label='Read Counts')
    plt.plot(smoothed_counts, ls='--', color='b')
    plt.axhline(filter_threshold, ls=':', color='r',
                label='Preliminary Filter\nThreshold = {}'.format(filter_threshold))
    plt.xlabel('UMIs')
    plt.ylabel('CCS Counts')
    plt.legend()
    plt.savefig(os.path.join(output_folder, f'{sample_name}_counts.png'), format='png', dpi=200)

    print("\tPreliminary count filter threshold: {}".format(filter_threshold))

    with open(os.path.join(output_folder, f'{sample_name}_umi_seq.txt'), 'w') as f_out:
        for count, umi in sorted_umis:
            if count >= filter_threshold:
                f_out.write(f'{umi}\n')


if __name__ == '__main__':
    fp_input = sys.argv[1]
    rt_primer_rc = sys.argv[2]
    output_dir = sys.argv[3]
    umi_len = int(sys.argv[4])
    main(fp_input, rt_primer_rc, output_dir, umi_len)
