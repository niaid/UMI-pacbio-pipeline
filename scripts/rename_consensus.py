from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys


if __name__ == "__main__":
    fp_in = sys.argv[1]
    sample_name = sys.argv[2]
    umi = sys.argv[3]
    count = sys.argv[4]

    with open(fp_in, 'r') as f:
        record = next(SimpleFastaParser(f))

    with open(fp_in, 'w') as f:
        f.write(f'>{sample_name}_{umi}_{count}\n')
        f.write(f'{record[1]}\n')