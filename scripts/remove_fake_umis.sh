#!/bin/bash

folder=$1
file=$2
n_cores=$3

out_after_cur=$folder"sgs/"
out_dir=$folder"sgs-unaligned/"
mkdir -p $out_dir
mkdir -p $out_after_cur

# Does fake umi removal
mkdir -p $folder"fake-umi-curation"
python3 $SCRIPT_DIR/umi_dedup.py $folder $file $n_cores
python3 $SCRIPT_DIR/inflection_removal.py $folder $file "fake-umi-curation" $infl_yes
cp $folder"fake-umi-curation/"$file"/"$file"_final.fasta" $out_after_cur
python3 $SCRIPT_DIR/down_select_post_cur.py $folder"sgs-prelim/"$file".fasta" $out_after_cur$file"_final.fasta" $out_dir$file"_unaligned.fasta"
