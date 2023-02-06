#!/bin/bash

RTprimer_RC=$1
project=$2
file=$3
barcode=$4

base=$project"sequences/"$file"/"
trim_inputFile=$base$barcode"_seq.fasta"
trimmed_outFile=$base"trimmed/"$barcode"_trimmed.fasta"

echo "Trimming RT primers..."
cutadapt -a $RTprimer_RC -o $trimmed_outFile $trim_inputFile > $project"sequences/"$file"/"usearch-cluster-log.txt

#=================================================================#

centroidFile=$base"centroid_usearch/"$barcode"_centroid.fasta"
consensusFile=$base"consensus_usearch/"$barcode"_consensus.fasta"
clusterFile=$base"cluster_usearch/"$barcode"_cluster.uc"

echo "Clustering..."
usearch -cluster_fast $trimmed_outFile -id 0.99 -centroids $centroidFile -consout $consensusFile -uc $clusterFile -threads 1 > $project"sequences/"$file"/"usearch-cluster-log.txt
