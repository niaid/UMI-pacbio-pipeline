"""
This script takes only sequences from final_ccs_reads that went through the fake umi process
"""

import sys
from Bio import SeqIO
from umi_dedup import strip_umi

def main(final_reads_path, post_curation_path, output_file):
	"""
	see module doc string

	Parameters:
	Fill in later
	"""

	# reads in the post curation file 
	post_curation_reads = SeqIO.parse(post_curation_path, "fasta")
	post_curation_dict = {strip_umi(x.id).split("_")[0]:x.id for x in post_curation_reads}

	final_reads = SeqIO.parse(final_reads_path, 'fasta')
	keeps = []
	for seq in final_reads:
		if any(item in seq.id for item in post_curation_dict.keys()):
			keeps.append(seq)

	for sequence in keeps:
		umi = strip_umi(sequence.id).split("_")[0]
		sequence.id = post_curation_dict[umi]
		sequence.name = post_curation_dict[umi]
		sequence.description = post_curation_dict[umi]

	with open(output_file, 'w') as handle:
		SeqIO.write(keeps, handle, "fasta")

if __name__ == '__main__':
	FINAL_READS = sys.argv[1]
	POST_CUR_PATH = sys.argv[2]
	OUTPUT_FILE = sys.argv[3]
	main(FINAL_READS, POST_CUR_PATH, OUTPUT_FILE)
