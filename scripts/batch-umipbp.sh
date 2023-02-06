#!/bin/bash

project=$1
config_file=$2
gen_dist=$3

# echo "Project: " $project > $project/log.txt
# echo "Configuration file: " $config_file >> $project/log.txt
# echo "Genetic distance: " $gen_dist >> $project/log.txt

# module unload all
# module load usearch
# module load "BLAST+"
# module load mafft
# module load Anaconda3/2020.07

# # Activate Conda environment with cutadapt, Biopython, joblib, etc. (see requirements.txt) - do not use LOCUS system cutadapt
# source activate pbpenv

ls $project | grep "fasta\$" | rev | cut -d "." -f2- | rev | while read -r file;
do
	echo $file
	./scripts/pacbio-pipeline-with-blast.sh $project $file $config_file $gen_dist
done
