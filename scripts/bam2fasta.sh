#!/bin/sh

# Job Name
#$ -N bam2fasta

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
# qsub bam2fasta.sh FILE_PATH OUTPUT
# OUTPUT must end in a "/"

FILE_PATH=$1
OUTPUT=$2


#/hpcdata/vrc/vrc1_data/smrtlink/jobs_root/cromwell-executions/pb_demux_ccs/7bbb9b83-00f9-4a7b-b68a-107da23b0396/call-demultiplex_barcodes/demultiplex_barcodes/6cd4cff9-10d7-4d69-8e54-0630d7592b88/call-lima/execution/
mkdir -p $OUTPUT

ls $FILE_PATH*.bam | rev | cut -d "." -f2- | cut -d "/" -f 1 | rev | while read -r file_in;
do
	/nethome/smrtanalysis2/smrtlink/smrtcmds/bin/bam2fasta -o $OUTPUT$file_in  $FILE_PATH$file_in.bam
done 

echo "Done!"

find $OUTPUT -size 0 -delete
rm $OUTPUT*remove*

#cat lima.lbc10--lbc10.consensusreadset.xml | grep -o -P  '(?<=ResourceId=").*(?=lima.lbc)' | uniq

