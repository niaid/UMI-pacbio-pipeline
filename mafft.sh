#!/bin/sh

# Job Name
#$ -N mafft

# Execute the script from the current working directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

#$ -S /bin/bash

# Send the output of the script to a directory called 'UGE-output' under current working directory (cwd)
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

# project: directory of fasta file, there should be / at the end.
# file: fasta file without format


project=$1
file=$2
#num_cor=$3

inputFile=$project$file".fasta"
outputFile=$project"mafft/"$file"_mafft.fasta"
refFile="reference-files/CoV_Reference_SEM_protein.fa"
module load mafft
mkdir -p $project"mafft"

#Progressive methods
# FFT-NS-1: * commonly used
#mafft --auto --retree 1 --maxiterate 0 $inputFile > $outputFile

#mafft --auto --add $inputFile $refFile > $outputFile

## extremely large sequences >200,000
#mafft --parttree --retree 1 $inputFile > $outputFile

#mafft --globalpair --thread 16 --maxiterate 1000 $inputFile > $outputFile

#########################################################################

## not good for large number of sequences good for around 5000

# FFT-NS-2: Improved version of the FFT-NS-1
mafft --auto --retree 2 $inputFile > $outputFile

# Iterative refinement method
#mafft --maxiterate 1000 $inputFile > $outputFile
