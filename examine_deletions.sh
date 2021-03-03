#!/bin/bash
# This script examines all the deletions in each sample

# Reads all bam files in the project directory
project=$1
outdir=$project"deletions/"

mkdir -p $outdir

module unload python
module load Anaconda3/2020.07
source activate umi-error

ls $project | grep "bam$\|bam\*$" | rev | cut -d "." -f 2- | rev | while read -r file_in;
do
 echo $project$file_in
 python examine_deletions.py "$project""$file_in".bam "$outdir$file_in.csv"
done

