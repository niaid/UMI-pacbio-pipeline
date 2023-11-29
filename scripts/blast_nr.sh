#!/bin/bash

project=$1
file=$2
db=$3
njobs=$4

inputfile=$project$file".fasta"
outputfile=$project"blastn-out/"$file".out"

# If you are querying a local db, point to the file inside the folder e.g.
# -db "/absolute/path/to/file.fa"
blastn -query $inputfile -db $db -num_threads $njobs -outfmt "7 qseqid evalue bitscore stitle" -out $outputfile