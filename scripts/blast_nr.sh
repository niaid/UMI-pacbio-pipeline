


#!/bin/sh

# Job Name
#$ -N blastn-nr

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
#$ -l h_vmem=100G

# qsub blast_nr.sh project file
# project: the directory where fasta file that blast removed is located, there should be / at the end.
# file: the name of fasta file without extenstion 

module unload all
module load "BLAST+"


project=$1
file=$2
inputfile=$project$file".fasta"
outputfile=$project"Blast/"$file"_out.out"

mkdir -p $project"Blast"

echo $inputfile

#blastn -db nt -query $inputfile -outfmt "6 stitle std qlen" -max_target_seqs 5 -out results.out -remote   
#blastn -db nt -query $inputfile -outfmt "6 qseqid evalue bitscore stitle" -out $outputfile -remote   
#blastn -db nt -query $inputfile -outfmt "7 qseqid evalue bitscore stitle" -max_target_seqs 1 -max_hsps 1 -out $outputfile -remote   
blastn -query $inputfile -db "/hpcdata/bio_data/blast_db_29SEP2020/nt" -num_threads 8 -outfmt "7 qseqid evalue bitscore stitle" -out $outputfile
