#!/bin/sh

# Job Name
#$ -N remove_fake_umis_nogen

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

# Tell the job your memory requirements
#$ -l h_vmem=5G

# Usage
# bash file
# 
project=$1
file=$2
n_cores=$3
infl_yes=$4

out_after_cur=$project"final_post_curation/nongen/"
mkdir -p $out_after_cur

# Does fake umi removal
module unload python
module load Anaconda3/2020.07
source activate umi-error
mkdir -p $project"fake-umi-curation-nogen"
python "umi_dedup_nogen.py" $project $file $n_cores
python "inflection_removal.py" $project $file "fake-umi-curation-nogen" $infl_yes
cp $project"fake-umi-curation-nogen/"$file"/"$file"_final.fasta" $out_after_cur
