import sys
from Bio import SeqIO


def main(fpf, fpu, output):

    # Read sequences
    seqs = dict()
    with open(fpf, 'r') as f:
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            seqs[record.name] = record.seq

    # Read UMIs
    umis = []
    with open(fpu, 'r') as f:
        for line in f:
            umi = line.strip()
            if umi:
                umis.append(umi)

    # Write sequences for each UMI to file
    for umi in umis:
        fp_out = f'{output}/{umi}.fasta'
        with open(fp_out, 'w') as f:
            for seq in seqs:
                if umi in seq:
                    f.write(f'>{seq}\n')
                    f.write(f'{seqs[seq]}\n')


if __name__ == '__main__':
    fp_fasta = sys.argv[1]
    fp_umis = sys.argv[2]
    output_folder = sys.argv[3]
    main(fp_fasta, fp_umis, output_folder)
