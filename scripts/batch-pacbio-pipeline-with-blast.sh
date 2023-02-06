#!/bin/sh

#$ -N pbptdev
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -o /hpcdata/vrc_vpds/radeckipe/dev/umipbp-dev/job-outputs/
#$ -m be
#$ -pe threaded 16
#$ -M radeckipe@nih.gov

# ./batch-pacbio-pipeline-with-blast.sh project config_file
# project: the folder where all the fasta files are. There must be a / at the end.
# config_file: This is the config script that loads all the variables.
#		Example config files are located in the folder config-scripts.
#		This allows change of primers and read length without editing code.
# This only operates on fasta files, so if there are other files/folders in the project directory, those are unchanged.

project=$1
config_file=$2

echo "Project: " $project > $project/project-log.txt
echo "Configuration file: " $config_file >> $project/project-log.txt
echo "Genetic distance: " $gen_dist >> $project/project-log.txt

module purge
module load uge
module load usearch
module load "BLAST+"
module load mafft
module load Anaconda3/2020.07
module load minimap2
module load samtools/1.13-GCC-4.8.4
module load bcftools/1.13-GCC-4.8.4
module load seqkit
export PATH=$PATH:/hpcdata/vrc_vpds/scripts/vsearch/vsearch-2.21.1-linux-x86_64-static/bin


ls $project | grep "fasta\$" | rev | cut -d "." -f2- | rev | while read -r file;
do
	echo $file
	module refresh
	. /hpcdata/vrc_vpds/scripts/umipbp-dev/scripts/pacbio-pipeline-with-blast.sh $project $file $config_file
done
