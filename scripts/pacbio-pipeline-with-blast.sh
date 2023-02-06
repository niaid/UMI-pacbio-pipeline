#!/bin/bash

project=$1
sample=$2
config_file=$3

input_file=$project$sample".fasta"
output_dir=$project"sgs/"

# -------------------------
# ----- Preliminaries -----
# -------------------------

# Recover script directory to ensure paths are resolved
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=$SCRIPT_DIR:$PATH
export PYTHONPATH=$SCRIPT_DIR:$PYTHONPATH
# Get specific version of pipeline currently being run, included git commit tag
version="$( get_version.py )-$( cd $SCRIPT_DIR &> /dev/null && git rev-parse --short HEAD )"

# Initialize QC log file
mkdir -p $project/qc
QCLOG=$project/qc/$sample-qc.txt
echo $(date) "umipbp version: " $version > $QCLOG
echo $(date) "PATH =" $PATH >> $QCLOG

# Load configuration details for sample
source $config_file
export PYTHONDONTWRITEBYTECODE=1

# Determine number of threads to use: if JOBS defined in config file, use JOBS number of threads. Otherwise,
# use as many as hardware allows, up to 16 maximum threads
nthreads=$(getconf _NPROCESSORS_ONLN)
echo $(date) "Number of available threads:" $nthreads >> $QCLOG
if [ -z ${JOBS+x} ]
then
  # If JOBS is not specified:
  njobs=$((nthreads < 16 ? nthreads : 16)) # Use up to 16 threads, if available
  echo $(date) "Maximum number of jobs to use:" $njobs >> $QCLOG
else
  # If JOBS -is- specified:
  njobs=$JOBS # Use manually-set number of threads
  echo $(date) "Maximum number of jobs to use:" $njobs >> $QCLOG

fi

echo $(date) "Starting to analyze sample:" $input_file >> $QCLOG
echo $(date) "Number of CCS reads:" $(grep -c ">" $input_file) >> $QCLOG

# Prepare output directory
mkdir -p $output_dir



# ----------------------------
# ----- Read Orientation -----
# ----------------------------

oriented_output_folder=$project"oriented"
mkdir -p $oriented_output_folder
oriented_file=$oriented_output_folder"/oriented_"$sample".fasta"

echo $(date) "Orienting reads with vsearch --orient" >> $QCLOG
#usearch -orient $input_file -db $referenceFile -fastaout $oriented_file -threads $njobs
vsearch --orient $input_file --db $referenceFile --fastaout $oriented_file --threads $njobs
echo $(date) "Reads oriented!" >> $QCLOG
echo $(date) "Number of CCS reads oriented:" $(grep -c ">" $oriented_file) >> $QCLOG



# -------------------------------
# ----- PCR Primer Trimming -----
# -------------------------------

# Trim PCR primers; sequences between $min_length and $max_length (after PCR primers removed) are kept
mkdir -p $project"/trimmed"

orientedFile=$oriented_output_folder"/oriented_"$sample".fasta"
untrimmedFile=$project"/trimmed/"$sample"_untrimmed.fasta"
trimmedFile=$project"/trimmed/"$sample"_trimmed.fasta"
shortFile=$project"/trimmed/"$sample"_short.fasta"
longFile=$project"/trimmed/"$sample"_long.fasta"

untrimmedForwardFile=$project"/trimmed/"$sample"_untrimmedF.fasta"
trimmedForwardFile=$project"/trimmed/"$sample"_trimmedF.fasta"

echo $(date) "Trimming primers and filtering by length with cutadapt" >> $QCLOG
cutadapt -j $njobs -e $error -g $forward --untrimmed-output $untrimmedForwardFile $orientedFile > $trimmedForwardFile
cutadapt -j $njobs -e $error -a $reversePrimer_RC --minimum-length $minLength --too-short-output $shortFile --maximum-length $maxLength --too-long-output $longFile --untrimmed-output $untrimmedFile $trimmedForwardFile > $trimmedFile
echo $(date) "Reads trimmed!" >> $QCLOG
echo $(date) "Number of CCS reads after trimming:" $(grep -c ">" $trimmedFile) >> $QCLOG


# -----------------------
# ----- UMI Binning -----
# -----------------------

mkdir -p $project"/umi-stats/"
echo "Using python executable at" $(which python)
echo $(date) "Binning by UMI..." >> $QCLOG
python3 $SCRIPT_DIR/bin_umis.py $trimmedFile $RTprimer_RC $project"/umi-stats/" $umiLength >> $QCLOG
echo $(date) "Done." >> $QCLOG
echo $(date) "Number of preliminary unique UMIs:" $(wc -l < $project"/umi-stats/"$sample"_counts.txt") >> $QCLOG

echo $(date) "Preliminary UMI stats:" >> $QCLOG
python3 $SCRIPT_DIR/umi_stats.py < $project"/umi-stats/"$sample"_counts.txt" >> $QCLOG


# This will make the folders required for the next phase. We separate this from the python script for parallel execution
rm -r $project"bins/"$sample
mkdir -p $project"bins/"$sample
mkdir -p $project"bins/"$sample"/with_rtp"
mkdir -p $project"bins/"$sample"/trimmed"
mkdir -p $project"bins/"$sample"/aligned"
mkdir -p $project"bins/"$sample"/vcf"
mkdir -p $project"bins/"$sample"/logs"
mkdir -p $project"bins/"$sample"/consensus"
mkdir -p $project"bins/"$sample"/clusters"

# This extracts umi reads into their own samples executed in parallel
echo $(date) "Extracting sequences from UMI bins..." >> $QCLOG
umi_seq=$project"umi-stats/"$sample"_umi_seq.txt"
fasta_umi_labeled=$project"umi-stats/"$sample"_UMI.fasta"
python3 $SCRIPT_DIR/extract_seqs.py $fasta_umi_labeled $umi_seq $project"bins/"$sample"/with_rtp"
echo $(date) "Done." >> $QCLOG


# -----------------------------------------------
# ----- Preliminary Consensus Determination -----
# -----------------------------------------------


echo $(date) "Determining UMI consensus sequences..." >> $QCLOG
xargs -n1 -P $(( njobs )) bash "flatten_bin.sh" $RTprimer_RC $project $sample $SCRIPT_DIR $referenceFile < $project"umi-stats/"$sample"_umi_seq.txt"
echo $(date) "Done." >> $QCLOG

#echo -e "\tNumber of bins with robust consensus:" $(cat $project"bins/"$sample"/logs/"*.log | grep -c "Robust") >> $QCLOG
#echo -e "\tNumber of bins with putative collision:" $(cat $project"bins/"$sample"/logs/"*.log | grep -c "Putative collisions") >> $QCLOG


# -----------------
# ----- BLAST -----
# -----------------

# Moves all the UMI consensus sequences to the same file, then BLASTs against reference database
consensus_folder=$project"bins/"$sample"/consensus"
blast_dir=$project"blast/"
blastn_dir=$blast_dir"blastn-out/"
matches_dir=$blast_dir"matches/"
prelim_sgs_dir=$project"sgs-prelim/"
mkdir -p $blast_dir
mkdir -p $blastn_dir
mkdir -p $matches_dir
mkdir -p $prelim_sgs_dir

# Copy consensus sequences to file
blast_input_file=$blast_dir"/"$sample".fasta"
echo "" > $blast_input_file

echo $(date) "Copying clusters..." >> $QCLOG
ls $consensus_folder | grep "fasta\$" | while read -r file;
do
	 cat $consensus_folder"/"$file >> $blast_input_file
done
echo $(date) "Done." >> $QCLOG

# Perform BLAST search and down-select matches
echo $(date) "BLASTing consensus sequences..." >> $QCLOG
bash "blast_nr.sh" $blast_dir $sample $dbFile $njobs
grep -i $term $blastn_dir$sample".out" | cut -f 1 | sort | uniq > $matches_dir$sample"_matches.txt"
echo $(date) "Down-selecting matches..." >> $QCLOG
python3 $SCRIPT_DIR/down_select_seqs.py $blast_dir$sample".fasta" $matches_dir$sample"_matches.txt" $prelim_sgs_dir$sample".fasta"
echo $(date) "Done." >> $QCLOG
echo $(date) "Preliminary number of SGSs:" $(grep -c ">" $prelim_sgs_dir$sample".fasta") >> $QCLOG

# ----------------------------
# ----- Fake UMI Removal -----
# ----------------------------

echo $(date) "Aligning preliminary SGSs..." >> $QCLOG
. mafft.sh $prelim_sgs_dir $sample $njobs
echo $(date) "Done." >> $QCLOG


mkdir -p $prelim_sgs_dir"sgs-consensus/"
mkdir -p $prelim_sgs_dir"fingerprints/"
mkdir -p $prelim_sgs_dir"fake-umi/"

# Identify PCR fingerprints + distinctive mutations
echo $(date) "Determining PCR fingerprints..." >> $QCLOG
python3 $SCRIPT_DIR/scan_fingerprints.py $project $sample $njobs
echo $(date) "Skipped." >> $QCLOG

# Do fake UMI removal
echo $(date) "Removing fake UMIs..." >> $QCLOG
python3 $SCRIPT_DIR/remove_fake_umis.py $project $sample > $prelim_sgs_dir"fake-umi/"$sample"_pairs.txt"
echo $(date) "Final number of SGSs:" $(grep -c ">" $prelim_sgs_dir"fake-umi/"$sample".fasta") >> $QCLOG


# ----------------------------------
# ----- Consensus Finalization -----
# ----------------------------------

pcr_dir=$prelim_sgs_dir"pcr-reversions/"
pcr_corrections=$pcr_dir"/"$sample"_pcr.txt"
mkdir -p $pcr_dir
touch $pcr_corrections

# PCR error removal: revert PCR errors that coincide with consensus sites to the consensus base
echo $(date) "Correcting PCR errors..." >> $QCLOG
python3 $SCRIPT_DIR/revert_pcr_errors.py $project $sample $output_dir$sample"_final.fasta" > $pcr_corrections
echo $(date) "Number of PCR errors corrected:" $(wc -l < $pcr_corrections) >> $QCLOG


# -----------------------------------------
# ----- Generate Final SGS Alignments -----
# -----------------------------------------

mkdir -p $output_dir"unaligned/"
mkdir -p $output_dir"reference-aligned/"

# Within-group alignment
echo $(date) "Aligning SGS against each other..." >> $QCLOG
mafft --thread $njobs $output_dir$sample"_final.fasta" > $output_dir$sample".fasta"
mv $output_dir$sample"_final.fasta" $output_dir"unaligned/"$sample".fasta"
echo $(date) "Done." >> $QCLOG

# Align against reference sequence
echo $(date) "Aligning SGS against reference..." >> $QCLOG
cat $referenceFile > $output_dir"reference-aligned/"$sample".fasta"
echo $'\n' >> $output_dir"reference-aligned/"$sample".fasta"
cat $output_dir$sample".fasta" >> $output_dir"reference-aligned/"$sample".fasta"
mafft --thread $njobs $output_dir"reference-aligned/"$sample".fasta" > $output_dir"reference-aligned/"$sample"_aligned.fasta"
rm $output_dir"reference-aligned/"$sample".fasta"
echo $(date) "Done." >> $QCLOG

echo $(date) "Finished processing "$sample"." >> $QCLOG
