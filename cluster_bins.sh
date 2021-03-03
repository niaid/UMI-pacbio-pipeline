#!/bin/sh
#$ -N trim
#$ -S /bin/bash
#$ -l quick,h_vmem=10G

#module unload all
module load cutadapt
module unload python

RTprimer_RC=$1
project=$2
file=$3
barcode=$4

base=$project"sequences/"$file"/"
trim_inputFile=$base$barcode"_seq.fasta"
trimmed_outFile=$base"trimmed/"$barcode"_trimmed.fasta"

cutadapt -a $RTprimer_RC -o $trimmed_outFile $trim_inputFile 

#=================================================================#

module load usearch
centroidFile=$base"centroid_usearch/"$barcode"_centroid.fasta"
consensusFile=$base"consensus_usearch/"$barcode"_consensus.fasta"
clusterFile=$base"cluster_usearch/"$barcode"_cluster.uc"

usearch -cluster_fast $trimmed_outFile -id 0.99 -centroids $centroidFile -consout $consensusFile -uc $clusterFile
