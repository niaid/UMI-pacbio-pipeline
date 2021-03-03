#!/bin/sh

# Job Name
#$ -N pacbio_pipeline-with-blast

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
#$ -l h_vmem=2G

# execute this script from the directory containing it.

# Run1: qsub -pe round 8 pacbio-pipeline-with-blast.sh project file config_file nogen  
#or
# Run2: qsub -pe round 8 pacbio-pipeline-with-blast.sh project file config_file   

# Parameters
# project: path to the folder containing the data file. Must end in a /. This must be absolute.
# file: This is the name of the fasta file to be processed. It should not contain the extension. It should also be inside the project folder
# config_file: This is the config script that loads all the variables. Must be an absolute path

module unload all

project=$1
file=$2
config_file=$3
gen_dist=$4
inputFile=$project$file".fasta"
oriented_outFile=$project"oriented"
output_dir=$project"final_ccs_reads"
tobe_blasted=$project"Tobe_blasted"

# Loads config file with all the variables
source $config_file

mkdir -p $output_dir


module load usearch
## orient the files with reference to hxb2 ##
mkdir -p $oriented_outFile
cd $oriented_outFile
echo "Correct orientation!"
usearch -orient $inputFile -db $referenceFile -fastaout $oriented_outFile"/orient_"$file".fasta" -notmatched $oriented_outFile"/noMatch_"$file".fasta"
cat $oriented_outFile"/orient_"$file".fasta" $oriented_outFile"/noMatch_"$file".fasta" > $oriented_outFile"/oriented_"$file".fasta"
rm $oriented_outFile"/noMatch_"$file".fasta" $oriented_outFile"/orient_"$file".fasta"
echo "Reads oriented!"
## -tabbedout text file can be added to get orientation of input sequence##

#=========================================================================#

module load cutadapt

#sequences shorter than 90% of the length of the insert are trimmed#
mkdir -p $project"/trimmed"

orientedFile=$oriented_outFile"/oriented_"$file".fasta"
untrimmedFile=$project"/trimmed/"$file"_untrimmed.fasta"
trimmedFile=$project"/trimmed/"$file"_trimmed.fasta"
shortFile=$project"/trimmed/"$file"_short.fasta"
longFile=$project"/trimmed/"$file"_long.fasta"

echo "Trim primers and filter by length"
cutadapt -e $error -g $forward -g $forwardPrimer_RC $orientedFile | cutadapt -e $error -a $reverse -a $reversePrimer_RC --minimum-length $minLength --too-short-output $shortFile --maximum-length $maxLength --too-long-output $longFile --untrimmed-output $untrimmedFile - > $trimmedFile

#=========================================================================#

 module unload python
 module load Anaconda3/2020.07
 source activate umi-error

 # Note the path of the python pacbio-pipeline. This shell command is expected to be run in the directory of the shell script.

 python "scripts/bin_by_UMI.py" $RTprimer_RC $project $file

 # This will make the folders required for the next phase. We separate this from the python script for parallel execution
 mkdir -p $project"sequences/"$file
 mkdir -p $project"sequences/"$file"/trimmed"
 mkdir -p $project"sequences/"$file"/centroid_usearch"
 mkdir -p $project"sequences/"$file"/consensus_usearch"
 mkdir -p $project"sequences/"$file"/cluster_usearch"

 # This extracts umi reads into their own files executed in parallel
 xargs -l -P 8 python "scripts/extract_seq.py" $project $file< $project"umi_stats/"$file"_umi_seq.txt"
 xargs -l -P 8 bash "scripts/cluster_bins.sh" $RTprimer_RC $project $file< $project"umi_stats/"$file"_umi_seq.txt"
 xargs -l python "scripts/read_uc_file.py" $project $file< $project"umi_stats/"$file"_umi_seq.txt"
 # next steps move reads into individual clusters for umis that could have collision
 # Does nothing for umis without collision
 cut -d "," -f 1 $project"umi_collision/"$file"collision.txt" > $project"umi_collision_all.tmp"
 if [ -s $project"umi_collision_all.tmp" ]
 then xargs -l python "scripts/umi_cluster_reads.py" $project $file< $project"umi_collision_all.tmp"
 fi
 # Cleans up useless files
 rm $project"umi_collision_all.tmp"
 find $project"error_insert" -size 0 -delete
 find $project"umi_collision" -size 0 -delete

 # Moves all the final reads to the same file
 mkdir -p $tobe_blasted
 cp $project"sequences/"$file"/cluster_with_read_counts/"$file"read.fasta" $tobe_blasted"/"$file"read.fasta"


 bash "scripts/blast_nr.sh" $project"Tobe_blasted/" $file"read"
 
 mkdir -p $project"Tobe_blasted/cov_headers"
 grep -i "coronavirus" $project"Tobe_blasted/Blast/"$file"read_out.out" | cut -f 1 | sort |uniq > $project"Tobe_blasted/cov_headers/"$file"_cov_headers.txt"
 python "scripts/down_select_seqs.py" $project"Tobe_blasted/"$file"read.fasta" $project"Tobe_blasted/cov_headers/"$file"_cov_headers.txt" $project"final_ccs_reads/"$file"read.fasta"  

 # Does fake umi removal
 uniq_umis=$(wc -l < $project"umi_stats/"$file"_umi_seq.txt")

 if [[ uniq_umis -le 500 ]]; then
 	num_cores=1
 elif [[ uniq_umis -gt 500 && uniq_umis -le 1500 ]]; then
 	num_cores=16
 elif [[ uniq_umis -gt 1500 && uniq_umis -le 2500 ]]; then
 	num_cores=32
 elif [[ uniq_umis -gt 2500 ]]; then
 	num_cores=64
 fi

 if [[ "$gen_dist" == "nogen" ]]
 then
 	qsub -pe round $num_cores "scripts/remove_fake_umis_nogen.sh" $project $file $num_cores
 else
 	#Do align and gen dist stuff
 	bash "scripts/mafft.sh" $project"final_ccs_reads/" $file"read"
 	qsub -pe round "$num_cores" "scripts/remove_fake_umis.sh" "$project" "$file" "$num_cores" 
 fi
 
