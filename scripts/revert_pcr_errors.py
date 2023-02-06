import sys
import os
from seqlib import gen_coord_map, get_counts, determine_consensus
from Bio.SeqIO.FastaIO import SimpleFastaParser



if __name__ == "__main__":

    project = sys.argv[1]
    sample = sys.argv[2]
    output = sys.argv[3]

    with open(os.path.join(project, f'sgs-prelim/fake-umi/{sample}.fasta'), 'r') as f:
        sgs_seqs = {record[0]: str(record[1]) for record in SimpleFastaParser(f)}

    with open(os.path.join(project, f'sgs-prelim/sgs-consensus/{sample}.fasta'), 'r') as f:
        sgs_cons = {record[0]: str(record[1]) for record in SimpleFastaParser(f)}
    sgs_cons = sgs_cons[f'{sample}_sgs-consensus']

    # final_seqs = dict()
    final_seqs = sgs_seqs
    with open(output, 'w') as f:
        for umi_record in sorted(final_seqs, key=lambda name: int(name.split('_')[-1]), reverse=True):
            f.write(f'>{umi_record}\n{final_seqs[umi_record]}\n')


    # for umi_record in sgs_seqs:
    #     umi = umi_record.split('_')[-2]
    #
    #     with open(os.path.join(project, f'bins/{sample}/aligned/{umi}.fasta'), 'r') as f:
    #         umi_seqs = {record[0]: str(record[1]) for record in SimpleFastaParser(f)}
    #
    #     counts = get_counts(umi_seqs)
    #     cons_dict = determine_consensus(counts, threshold=0)
    #     umi_cons = [cons_dict[i] for i in range(len(counts))]
    #
    #     lg_map = gen_coord_map(umi_cons, sgs_seqs[umi_record])
    #
    #     for i, counter in enumerate(counts):
    #         mc = counter.most_common(2)
    #         if mc[0][1] / sum(counter.values()) < 0.8:
    #             if lg_map[i] and sgs_cons[lg_map[i]] in (mc[0][0], mc[1][0]):
    #                 if umi_cons[i] != sgs_cons[lg_map[i]]:
    #                     print(umi, 'correction:', i, lg_map[i], umi_cons[i], sgs_cons[lg_map[i]])
    #                     umi_cons[i] = sgs_cons[lg_map[i]]
    #
    #     umi_cons = ''.join(umi_cons).replace('-', '')
    #     final_seqs[umi_record] = umi_cons
    #
    # with open(output, 'w') as f:
    #     for umi_record in sorted(final_seqs, key=lambda name: int(name.split('_')[-1]), reverse=True):
    #         f.write(f'>{umi_record}\n{final_seqs[umi_record]}\n')




