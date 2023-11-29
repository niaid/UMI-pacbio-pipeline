#!/bin/bash

folder=$1
name=$2
njobs=$3

inputFile=$folder$name".fasta"
outputFile=$folder"mafft/"$name".fasta"
mkdir -p $folder"mafft"

mafft --auto --thread $njobs --retree 2 $inputFile > $outputFile

