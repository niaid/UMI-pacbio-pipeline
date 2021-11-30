#!/bin/sh

# Job Name
#$ -N batch-process-mafft-multi-core

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


# to run execute the command 
# navigate to the folder containing the batch-mafft then run the following command
# qsub batch_mafft.sh project 
# project: path to a directory  with fasta files in it. There must be a / at the end.


project=$1


ls $project | grep "fasta\$" | rev | cut -d "." -f2- | rev | while read -r file;
do
	qsub scripts/mafft.sh $project $file
done
