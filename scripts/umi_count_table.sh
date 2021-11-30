#!/bin/sh

# This script counts the number of unique UMIs for each cell and make a csv file.
# to run execute the command 
# navigate to the folder containing the *_umi_seq
# project: path to a directory with .txt file. There must be a / at the end.


project=$1
file_in=$project$"*_umi_seq.txt"

echo 'sample name,unique umi count' > $project"umi_Stat.csv"

for f in $file_in
do
	sample=$(echo $f | rev | cut -d "/" -f 1 | rev | cut -d "_" -f 1)
	umiCount=$(wc -l < $f)
	echo $sample","$umiCount >> $project"umi_Stat.csv"
done 

echo "Done!"

