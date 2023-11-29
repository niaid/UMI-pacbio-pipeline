#!/bin/bash

RTprimer_RC=$1
project=$2
file=$3
SCRIPT_DIR=$4
rfile=$5
umi=$6


bin_folder=$project"bins/"$file"/"
raw_bin_file=$bin_folder"/with_rtp/"$umi".fasta"
bin_file=$bin_folder"/trimmed/"$umi".fasta"
bin_cluster_file=$bin_folder"/trimmed/"$umi"-cc.fasta"
aligned_bin_file=$bin_folder"/aligned/"$umi".sam"
aligned_sorted_bin_file=$bin_folder"/aligned/"$umi"-sorted.sam"
vcf_file=$bin_folder"/vcf/"$umi".vcf"
cluster_stats_file=$bin_folder"/clusters/"$umi".uc"
cluster_cons_file=$bin_folder"/clusters/"$umi"-clusters.fasta"
bin_cons_file=$bin_folder"/clusters/"$umi"-cons.fasta"
consensus_file=$bin_folder"/consensus/"$umi".fasta"

LOG_FLATTEN=$bin_folder"/logs/"$umi".log"

# ===== Trimming =====

echo $(date) "Trimming RT primers..." > $LOG_FLATTEN
cutadapt -a $RTprimer_RC -o $bin_file $raw_bin_file >> $LOG_FLATTEN 2> /dev/null


# ===== Clustering =====
echo $(date) "Clustering bin reads..." >> $LOG_FLATTEN
vsearch -cluster_fast $raw_bin_file -id 0.99 -iddef 1 -uc $cluster_stats_file -threads 1 -consout $cluster_cons_file >> $LOG_FLATTEN 2> /dev/null
python3 $SCRIPT_DIR/prepare_bin_clusters.py $cluster_stats_file $bin_file $cluster_cons_file $bin_cons_file $bin_cluster_file

if [ ! -e "$bin_cluster_file" ]; then
    echo "Clustering failed to generate a suitable bin for" $umi >> $LOG_FLATTEN
    exit
else
    echo "Clustering successful for" $umi "(Number of reads:"$(grep -c '>' $bin_cluster_file)")" >> $LOG_FLATTEN
fi


# ===== Aligning =====
echo $(date) "Assembling consensus..." >> $LOG_FLATTEN
minimap2 -t 1 -ax splice:hq $bin_cons_file $bin_cluster_file > $aligned_bin_file 2> /dev/null
samtools sort $aligned_bin_file > $aligned_sorted_bin_file 2> /dev/null
samtools index $aligned_sorted_bin_file 2> /dev/null
bcftools mpileup -d 10000 -X pacbio-ccs -f $bin_cons_file $aligned_sorted_bin_file 2> /dev/null | bcftools call --threads 1 -c --ploidy 1 > $vcf_file 2> /dev/null
bgzip $vcf_file >> $LOG_FLATTEN 2> /dev/null
bcftools index --threads 1 $vcf_file.gz >> $LOG_FLATTEN 2> /dev/null
bcftools consensus -a N -f $bin_cons_file -o $consensus_file-tmp $vcf_file.gz >> $LOG_FLATTEN 2> /dev/null
seqkit -is replace -p "n" -r "" $consensus_file-tmp > $consensus_file 2> /dev/null
rm -f $consensus_file-tmp
python3 $SCRIPT_DIR/rename_consensus.py $consensus_file $file $umi $(samtools view -c $aligned_bin_file)
echo "Finished processing bin for:" $umi >> $LOG_FLATTEN

# ===== Aligning =====
#
#echo $(date) "Aligning trimmed reads..."  >> $LOG_FLATTEN
#mafft --thread 2 $bin_file > $aligned_bin_file 2>> $LOG_FLATTEN
#
## ===== Consensus Determination =====
#
#echo $(date) "Flattening alignment..."  >> $LOG_FLATTEN
#python3 $SCRIPT_DIR/determine_consensus.py $file $aligned_bin_file $consensus_file >> $LOG_FLATTEN
