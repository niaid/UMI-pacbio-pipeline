#!/bin/sh

# Job Name
#$ -N batch-revert-mutations

# Execute the script from the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

#$ -S /bin/bash

# Send the output of the script to a directory called 'UGE-output' uder current working directory (cwd)
  if [ ! -d "UGE-output" ]; then #Create output directory in case it does NOT exist
      mkdir UGE-output
  fi
#$ -o UGE-output/

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M ellie.bayatmokhtari@nih.gov

# Tell the job your memory requirements
#$ -l h_vmem=5G

# Usage


INPUT_FOLDER=$1
VCF_FLAG=$2

OUTPUT=$INPUT_FOLDER"reverted_mutations/"
mkdir -p $OUTPUT

module load samtools
module unload python
module load Anaconda3/2020.07
source activate umi-error

if [[ "$VCF_FLAG" == "True" ]]
then
	ls $INPUT_FOLDER | grep "bam$" | cut -d "_" -f 1 | while read -r file_in;
	do
	    in_bam=$(ls $INPUT_FOLDER | grep $file_in | grep "bam" | grep -v "bai")
		VCF_FILE=$(ls $INPUT_FOLDER | grep $file_in | grep ".vcf")

		python "scripts/revert_mutations_bam.py" "$INPUT_FOLDER$in_bam" "$INPUT_FOLDER$VCF_FILE" $VCF_FLAG > $OUTPUT$file_in"_reverted.bam" 2> $OUTPUT$file_in"deleted_inserts.txt"

		samtools index $OUTPUT$file_in"_reverted.bam"
		samtools fasta $OUTPUT$file_in"_reverted.bam" > $OUTPUT$file_in"_reverted.fasta"
	done
else
	var_file=$(ls $INPUT_FOLDER | grep ".xlsx$" | grep -v '^~')

	python make_variant_dict.py "$INPUT_FOLDER$var_file" $OUTPUT"variant_dict.pickle"  
	ls $INPUT_FOLDER | grep "bam$" | cut -d "_" -f 1 | while read -r file_in;
	do
		in_bam=$(ls $INPUT_FOLDER | grep $file_in | grep "bam" | grep -v "bai")
		
		python "scripts/revert_mutations_bam.py" "$INPUT_FOLDER$in_bam" $OUTPUT"variant_dict.pickle" > $OUTPUT$file_in"_reverted.bam" 2> $OUTPUT$file_in"deleted_inserts.txt"

		samtools index $OUTPUT$file_in"_reverted.bam"
		samtools fasta $OUTPUT$file_in"_reverted.bam" > $OUTPUT$file_in"_reverted.fasta"
	done 
fi

