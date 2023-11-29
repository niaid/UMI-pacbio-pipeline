"""
revert-mutations-reference.py [file]

Reverts technical errors to consensus and computes haplotypes.

Inputs:
    Fasta file of aligned SGS, aligned to a reference sequence which is the first sequence

Outputs:
    SGS with called variants kept, variant below threshold reverted ([file]-reverted.fa)
    Haplotype sequences ([file]-haplotypes.fa)

Copyright (c) 2022-2023 National Institutes of Health
Written by Pierce Radecki

"""

import sys
import scipy.stats
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter
import re

RT_SNV_ERROR_RATE = 1e-4
MINIMUM_INDEL_COUNT = 3

regex_deletion = re.compile(r'-+')


def compute_threshold(n_sgs, ins_len, p, rt_err_rate=1e-4):
    pts = np.arange(0, n_sgs + 1)
    with np.errstate(divide='ignore'):
        cdf = 1 - scipy.stats.binom(n=n_sgs, p=rt_err_rate).cdf(pts - 1) ** ins_len
    thresh = np.where(cdf < p)[0][0] + 1
    return thresh


def determine_consensus(counts, threshold=0.7):
    cons = dict()
    for i, counter in enumerate(counts):
        c = counter.most_common(1)[0][1] / sum(counter.values())
        if c >= threshold:
            cons[i] = counter.most_common(1)[0][0]
    return cons


def parse_deletions(seq_dict, ref_len, minimum_indel_count=3):
    ns = len(seq_dict)
    deletions = {}
    for s in seq_dict.values():
        seq_deletions = []
        matches = regex_deletion.finditer(''.join(s))
        for match in matches:
            seq_deletions.append((match.span()[0], match.span()[1]))
        if seq_deletions:
            for d in seq_deletions:
                sd = str(d)
                if sd in deletions:
                    deletions[sd]['count'] += 1
                else:
                    deletions[sd] = {'tuple': d, 'count': 1}

    deletions_final = set()
    for d in deletions:
        dt = deletions[d]['tuple']
        # Deletions of 1 or 2 bases must be seen at >20% to be considered
        if dt[1] - dt[0] < 3 or dt[0] > ref_len - 4:
            if deletions[d]['count'] >= minimum_indel_count and deletions[d]['count'] / ns > 0.2:
                deletions_final.add(deletions[d]['tuple'])
        else:
            # Deletions >=3 bases must be seen more than the minimum indel count threshold
            if deletions[d]['count'] >= minimum_indel_count:
                deletions_final.add(deletions[d]['tuple'])
    return deletions_final


def revert_mutations(seq_dict, ref, minimum_snv_count, minimum_indel_count=3):
    ref_end = len(ref.rstrip('-'))
    ref_start = len(ref) - len(ref.lstrip('-'))
    ccs_len = len(next(iter(seq_dict.values())))

    for index in range(ccs_len):
        if index < ref_start or index > ref_end - 10:
            # Ignore SGS sequences outside of reference sequence; omit last 10 bases where technical errors are common
            for s in seq_dict:
                seq_dict[s][index] = ref[index]
            continue

    counts = [Counter(''.join((seq_dict[s][j] for s in seq_dict))) for j in range(ccs_len)]  # Base counts per position
    cons = determine_consensus(counts, threshold=0)
    for index, counter in enumerate(counts):
        ns = sum(counter.values())
        for char in counter:
            # Check for insertions
            if counter.most_common(1)[0][0] == '-' and char != '-' and counter[char] / sum(counter.values()) < 0.2:
                for s in seq_dict:
                    if seq_dict[s][index] == char and index in cons:
                        seq_dict[s][index] = cons[index]  # Revert base
                continue
            if char != '-' and counter[char] < minimum_snv_count:
                for s in seq_dict:
                    if seq_dict[s][index] == char and index in cons:
                        seq_dict[s][index] = cons[index]

    # Check and revert or keep deletions: use three passes
    deletions = parse_deletions(seq_dict, ref_end, minimum_indel_count=minimum_indel_count)

    for s in seq_dict:
        matches = regex_deletion.finditer(''.join(seq_dict[s]))
        for match in matches:
            if (match.span()[0], match.span()[1]) in deletions:
                continue  # Keep deletion
            else:
                for i in range(match.span()[0], match.span()[1]):
                    seq_dict[s][i] = cons[i]  # Revert

    # Second pass
    deletions_final = parse_deletions(seq_dict, ref_end, minimum_indel_count=minimum_indel_count)

    counts_no_gaps = counts.copy()
    for counter in counts_no_gaps:
        counter['-'] = 0.5  # Hack to avoid calling '-' as consensus unless it's the only base at position
    cons_no_gaps = determine_consensus(counts_no_gaps, threshold=0)

    for s in seq_dict:
        matches = regex_deletion.finditer(''.join(seq_dict[s]))
        for match in matches:
            if (match.span()[0], match.span()[1]) in deletions_final:
                continue  # Keep deletion
            else:
                for i in range(match.span()[0], match.span()[1]):
                    seq_dict[s][i] = cons_no_gaps[i]  # Revert

    counts_final = [Counter(''.join((seq_dict[s][j] for s in seq_dict))) for j in range(ccs_len)]
    empty_cols = []
    for i, counter in enumerate(counts_final):
        if counter['-'] == sum(counter.values()) and ref[i] == '-':
            empty_cols.append(i)

    # Get rid of any empty columns remaining
    seqs_final = dict()
    for s in seq_dict:
        seqs_final[s] = [a for i, a in enumerate(seq_dict[s]) if i not in empty_cols]
    ref = [a for i, a in enumerate(ref) if i not in empty_cols]

    seq_len = len(next(iter(seqs_final.values())))
    counts_final = [Counter(''.join((seqs_final[s][j] for s in seqs_final))) for j in range(seq_len)]
    cons_final = determine_consensus(counts_final, threshold=0)

    # Final pass of deletion examination
    deletions_final = parse_deletions(seqs_final, len(ref), minimum_indel_count=minimum_indel_count)

    for s in seqs_final:
        matches = regex_deletion.finditer(''.join(seqs_final[s]))
        for match in matches:
            if (match.span()[0], match.span()[1]) in deletions_final:
                continue  # Keep deletion
            else:
                for i in range(match.span()[0], match.span()[1]):
                    seqs_final[s][i] = cons_final[i]  # Revert

    # Generate polished final sequences, reference
    seqs_final_final = dict()
    for s in seqs_final:
        seqs_final_final[s] = ''.join([a for i, a in enumerate(seqs_final[s])])
    ref = ''.join([a for i, a in enumerate(ref)])

    return seqs_final_final, ref


def get_haplotypes(seqs_dict):
    haps = dict()
    for sl in seqs_dict:
        s = ''.join(seqs_dict[sl])
        if s in haps:
            haps[s] += 1
        else:
            haps[s] = 1
    haps = {s: haps[s] for s in haps if haps[s] >= 2}  # Minimum count of 2
    return haps


if __name__ == "__main__":

    fasta_in = sys.argv[1]  # Fasta file of SGS

    # Read FASTA file
    with open(fasta_in, 'r') as f:
        seqs = {record[0]: [a for a in record[1]] for record in SimpleFastaParser(f)}

    ref_name = str(next(iter(seqs.keys())))
    print('Reference sequence:', ref_name)
    ref_seq = ''.join(seqs[ref_name].copy())
    del seqs[ref_name]

    # If we have fewer than 5 sequence, don't try to call haplotypes
    if len(seqs) < 5:
        sys.exit()

    min_snv_count = compute_threshold(len(seqs), len(ref_seq), 0.01)  # FDR of 0.01 per sample
    print("Minimum SNV count: ", min_snv_count)
    seqs_reverted, ref_seq = revert_mutations(seqs, ref_seq, min_snv_count)  # Revert technical errors

    # Write reverted sequences to [file]-reverted.fa
    with open(fasta_in.replace('.fa', '-reverted.fa'), 'w') as f:
        f.write(f'>{ref_name}\n{ref_seq}\n')
        for seq in seqs_reverted:
            nucs = ''.join(seqs_reverted[seq])
            f.write(f'>{seq}\n{nucs}\n')

    # Get haplotypes (identify and count unique sequences)
    haplotypes = get_haplotypes(seqs_reverted)

    # Write haplotypes to [file]-haplotypes.fa
    with open(fasta_in.replace('.fa', '-haplotypes.fa'), 'w') as f:
        f.write(f'>{ref_name}\n{ref_seq}\n')
        for n, seq in enumerate(sorted(haplotypes, key=lambda s: haplotypes[s], reverse=True)):
            nucs = ''.join(seq)
            f.write(f'>H{n + 1}_{haplotypes[seq]}\n{nucs}\n')
