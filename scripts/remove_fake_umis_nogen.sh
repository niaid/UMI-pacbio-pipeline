#!/bin/bash

folder=$1
file=$2
n_cores=$3
infl_yes=$4

out_after_cur=$folder"sgs/"
mkdir -p $out_after_cur

# Does fake umi removal
mkdir -p $folder"fake-umi-curation"
python $SCRIPT_DIR/umi_dedup_nogen.py $folder $file $n_cores
python $SCRIPT_DIR/inflection_removal.py $folder $file "fake-umi-curation" $infl_yes
cp $folder"fake-umi-curation/"$file"/"$file"_final.fasta" $out_after_cur
