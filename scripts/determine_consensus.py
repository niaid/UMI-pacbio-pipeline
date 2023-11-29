
import numpy as np
import sys
from Bio import SeqIO
from collections import Counter



def count_seqs(seqs):
    return [Counter(''.join(seq[i] for seq in seqs)) for i in range(len(seqs[0]))]


def determine_consensus(fp_input, fp_output):

    with open(fp_input, 'r') as handle:
        seqs = [record.seq for record in SeqIO.parse(handle, 'fasta')]

    counts = count_seqs(seqs)  # Position-wise counts

    bin_length = len(counts)

    cons = ['-'] * bin_length  # Consensus sequence
    p_cons = np.zeros(bin_length)  # Proportion of most common base

    p_collisions = []  # List to hold putative collisions, early-cycle PCR errors, etc.

    for i, counter in enumerate(counts):
        cons[i] = counter.most_common()[0][0]
        p_cons[i] = counter.most_common()[0][1] / sum(counter.values())
        if p_cons[i] < 0.75:
            p_collisions.append((i+1, counter.most_common()))

    if p_collisions:
        print("Putative collisions detected:")
        print("Index (1-indexed)\tNucleotide counts")
        for pc in p_collisions:
            print(f'{pc[0]}\t{pc[1]}')
        print("End of putative collisions.")
    else:
        print("Robust consensus sequence identified.")

    cons = ''.join(cons).replace('-', '')

    return cons, p_cons, len(seqs)


if __name__ == '__main__':
    sample_name = sys.argv[1]
    fp_aligned = sys.argv[2]
    fp_out = sys.argv[3]

    umi = fp_aligned.split(r'/')[-1].split('.fa')[0]

    # Determine consensus sequence
    cs, pcs, bin_size = determine_consensus(fp_aligned, fp_out)

    with open(fp_out, 'w') as f:
        f.write(f'>{sample_name}_{umi}_{bin_size}\n')
        f.write(f'{cs}\n')

