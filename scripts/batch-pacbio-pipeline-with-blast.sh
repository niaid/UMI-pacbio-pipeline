#!/bin/sh

# Job Name
#$ -N batch-pipeline-with-blast

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

# Tell the job your memory requirements
#$ -l h_vmem=2G

# This batch processes all the fasta files in a single directory
# to run execute the command 
# navigate to the folder contraining the batch-script then run the following command
# qsub batch-pipeline.sh project config_file
# project: the folder where all the fasta files are. There must be a / at the end.
# config_file: This is the config script that loads all the variables.
#		Config files are located in the folder config-scipts.
#		This allows change of primers and read length without editing code.
# This only operates on fasta files, so if there are other files/folders in the project directory, those are unchanged.

project=$1
config_file=$2
gen_dist=$3

ls $project | grep "fasta\$" | rev | cut -d "." -f2- | rev | while read -r file;
do
	echo $file
	qsub -pe round 8 scripts/pacbio-pipeline-with-blast.sh $project $file $config_file $gen_dist
done
