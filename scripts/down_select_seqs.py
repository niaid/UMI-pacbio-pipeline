import os
import sys
from Bio import SeqIO

def main(fp_in, fp_seq_names, fp_out):
    fasta_file = list(SeqIO.parse(fp_in, "fasta"))
    seqs = set()
    with open(fp_seq_names) as f:
        for line in f:
            line = line.strip()
            if line:
                seqs.add(line)

    seqs2keep = []
    for seq in fasta_file:
        if seq.id in seqs:
            seqs2keep.append(seq)

    if os.path.exists(fp_out):
        os.remove(fp_out)
    with open(fp_out, "w") as f:
        SeqIO.write(seqs2keep, f, "fasta")



if __name__ == '__main__':
    fp_in = sys.argv[1]
    fp_seq_names = sys.argv[2]
    fp_out = sys.argv[3]
    main(fp_in, fp_seq_names, fp_out)

