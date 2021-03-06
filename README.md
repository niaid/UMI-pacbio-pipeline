# Pacbio UMI Analysis Pipeline
This pipeline was developed in Virus Persistence and Dynamics Section (VPDS) of the Immunology Laboratory (IL) in the Vaccine Research Center (VRC) at the National Institutes of Health (NIH).

## Requirements

This pipeline requires a python environment. A conda environment is recommended, and a requirements.txt is included listing all required modules.

### Note 
You will need to change the loading of your python environment, inside any scripts.
	   	 
### pacbio-pipeline-with-blast.sh

1. Orients the fasta files with reference to `*.fasta`, output generated to `oriented` folder.

2. Trim the forward  and reverse (PCR) primers sequence using cutadapt. This step also does an insert filtering. Inserts shorter than 90% of the length of the desired insert length are trimmed. The output is generated to `trimmed` folder.

3. Script `bin_by_UMI.py` gets called and generates results under the `umi_stats` folder. In this folder `*_umi_seq.txt` gives the list of umi sequences extracted from the data.

4. Script `extract_seq.py` gets called and extracts reads to `sequences` folder based on associated UMI sequence.

5. Script `cluster_bins.sh` gets called to get consensus sequence from each umi bin. The output is generated to `sequences` folder.

6. Script `read_uc_file.py` gets called to get a fasta file with consensus sequence with highest read abundance from each UMI bin. The output is generated in `sequences/cluster_with_read_counts` folder.  It then removes insert errors. The output of error reads is saved in `error_insert` folders. 

7. Script `umi_cluster_reads.py` gets called to extract and save insert error in `error_insert` folder. 

8. CCS reads obtained from step 1 to 7 are saved in `Tobe_blasted` folder for blasting against reference.

9. Script `blast_nr.sh` gets called to blast sequences in each fasta file against local blast database. The results are saved in two folders: 1. `Tobe_blasted/Blast/.out` which contains the blast results and `Tobe_blasted/cov_headers/.txt` which grab the headers with contains the word we are interested (for example. coronavirus). 

10. Script `down_select_seqs.py` gets called to subset the sequences with headers from `Tobe_blasted/cov_headers` folder for each fasta file. The output is saved in `final_ccs_reads` folder. 

11. Script `remove_fake_umis_nongen.sh` or `remove_fake_umis.sh` get called to detect and set aside artifactual UMIs. The output is saved in either `fake-umi-curation-nogen` or `fake-umi-curation`, respectively.

12. Final output is saved in  `final_post_curation/nongen` folder or `final_post_curation/gen` folder based on the choice in 11 as a `*.fasta` file.

## Usage

How to execute the pipeline:

1. Navigate to the pacbio-pipeline folder.

2. To run in batch execute the following command: 

Run1: 

	scripts/batch-pacbio-pipeline-with-blast.sh <path-to-raw-fasta-files-directory> <path-to-config-script>  nogen   

or

Run2: 

	scripts/batch-pacbio-pipeline-with-blast.sh <path-to-raw-fasta-files-directory> <path-to-config-script> config_file
 

where `path-to-raw-fasta-files-directory` is the folder where all the fasta files are located. There must be a `/` at the end. The path to the `path-to-raw-fasta-files-directory` folder needs to be an absolute path. `config_file` is the folder where the primers, error, length information are located.
`gen_dist` is nogen in case of excluding genetic distance, or nothing which include genetic distance in network-based filtering.

Note: This only operates on fasta files, so if there are other files/folders in the project directory, those are unchanged.

3. To run in individual fasta file execute the following command:

Run1:

 	 qsub -pe round 8 scripts/pacbio-pipeline-with-blast.sh project file config_file nogen

or

Run2: 

 	 qsub -pe round 8 scripts/pacbio-pipeline-with-blast.sh project file config_file  


### Example

	qsub scripts/batch-pacbio-pipeline-with-blast.sh <path-to-raw-fasta-files-directory> <path-to-config-script> nogen

or 

	qsub scripts/batch-pacbio-pipeline-with-blast.sh <path-to-raw-fasta-files-directory> <path-to-config-script>



# Revert Mutations by Erroneous to Consensus Pipeline:

This pipeline, conserves the positions of real variants detected from variant caller and revert the rest of the bases back to consensus. <!---This is now a pypi package see [rev-seqs](https://pypi.org/project/rev-seqs/)-->

## Components:

### make_variant_dict.py

This script reads variant excel file (obtained from variant caller) and turns this into a dictionary and save this dictionary in pickle format to be read later. The output contains all the real variants.

### revert_mutations_bam.py

This script reads either pickle dictionary or .vcf file and bam/bai file and returns bam file where each file was reconstructed by conserving real variants and revert the rest of the bases to consensus. 

### batch-revert_mutations.sh

This script reads .xlsx( if True .vcf), .bam, .bai files and save the output as .bam format in `reverted_mutations` folder. It then use samtools to extract index and fasta files for each .bam file. The results are saved in `reverted_mutations` folder. 

## Usage

How to execute:

1. Navigate to the pacbio-pipeline folder.

2. To run in batch execute the following command: 

Run1: 

	scripts/batch-revert_mutations.sh <path-to-bam-files-directory>  True  

or

Run2: 

	scripts/batch-revert_mutations.sh <path-to-bam-files-directory> 
 

where `path-to-bam-files-directory` is the folder where all the bam/.bai/.vcf files are located. There must be a `/` at the end. 

If `True`, it takes .vcf file for variant caller. 

If none, it take excel file for variants information in .xlsx format.


### Example

	qsub scripts/batch-revert_mutations.sh <path-to-bam-files-directory> True

or 

	qsub scripts/batch-revert_mutations.sh <path-to-bam-files-directory>

## Pipeline components:
	
## bin_by_UMI.py
This script has functions that read the trimmed files and looks for 8 base primerID sequence after RT primer. Each primerID (UMI) with associated ID and sequence is written to a `*.csv` file. UMIs and their frequency of occurances are written to a `*.txt`. UMI sequences with read counts above the inflection point in UMI count distribution are preserved and formatted file as : >seqID_UMI \n sequence is written to a `*.fasta` file. 

## extract_seq.py
This script extracts reads to `sequences` folder based on associated UMI sequence. It collects all the sequences with the same UMI and place them in one fasta file formated as UMI_seq.fasta

## cluster_bins.sh
This script first uses cutadapt to trim the RT primer and the 8 base umi sequence and then uses usearch -cluster_fast to generate consensus sequence based on 99% identity. This script generates a `trimmed` folder that has trimmed inserts left after removal of RT primer and 8 base UMI sequence. It also generates `centrioid_usearch` folder that has centroid sequences from usearch, `cluster_usearch` folder that has the `*.uc` files generated from running usearch. These `*.uc` files are used to get the read abundance for each umi sequence.

## read_uc_file.py
This script reads the `*.uc` files generated by running usearch and then calculates the read abundance associated with consensus sequence and writes to `cluster_stats` folder. The script then reads the consensus with highest read abundance and writes that to a folder `cluster_with_read_counts`. This folder has a final fasta file that can be used for multiple sequence alignment. It then removes insert error based on a sets of criteria. The output of error reads is saved in `error_insert` folders. CCS reads are saved in `final_ccs_reads` folder.

## umi_cluster_reads.py 
This script gets called to extract and save possible insert error in `error_insert` folders. 

## check_umi_counts.sh
This script checks if the number of unique UMIs in final output matches the number of UMIs in final reads + insert_error. If there is a mismatch it save the result as `*.txt` file in `project` folder, else it prints nothing.  

## remove_fake_umis_nogen.sh 
This script uses two functions: `umi_dedup_nogen.py`: this script locates and set aside fake UMIs considering edit distance 1, and a_n>2n_b-1 criteria. The output is fake UMIs that are saved in `fake-umi-curation-nogen` folder.`inflection_removal.py`: this script removes low count UMIs based on counts below inflection point or knee point. The choice of inflection/knee point is deterimed by True (use inflection point) and without (True) knee point. Final output is saved in  `final_post_curation/nongen` folder in `*.fasta` formatted files.

## remove_fake_umis.sh 
This script uses three functions: `umi_dedup.py`: this script locates and set aside fake UMIs considering genetic distance 0, edit distance 1, and a_n>2n_b-1 criteria.The output is fake UMIs that are saved in `fake-umi-curation` folder. `inflection_removal.py`: this script removes low count UMIs based on counts below inflection point or knee point. The choice of inflection/knee point is deterimed by True (use inflection point) and without (True) knee point. `down_select_post_cure.py`: this script calls final_ccs_reads prior to fake UMI curation and extract the headers and match with the headers in `final_post_curation/gen` files and select matched headers and saved unaligned final output in `final_unaligned_post_cur` folder. Final output is saved in  `final_post_curation/gen` folder in `*.fasta` formatted files.
	    	    
## blast_nr.sh
This script uses each fasta file from `Tobe_blasted` folder and blast it against local blast db. The output is printed to a file with .out format for each fasta file. This file contains information about each sequence and what it is blasted to. The output is saved in `Tobe_blasted/Blast` folder. Following this step, for each .out file, the `grep` tool is called to extract the sequences' headers that contains "coronavirus (in case of CoV2)" for each fasta file and save these headers as a `*.txt` file in `Tobe_blasted/cov_headers` folder. 

## down_select_seqs.py
This script uses header files from `Tobe_blasted/cov_headers` folder to subset sequences and write them as a fasta file in `file_ccs_reads` folder. 


## Citing UMI-pacbio-pipeline
If you use UMI-pacbio-pipeline in your work, please cite:

[Ko SH, Bayat Mokhtari E, Mudvari P, Stein S, Stringham CD, et al. (2021) High-throughput, single-copy sequencing reveals SARS-CoV-2 spike variants coincident with mounting humoral immunity during acute COVID-19. PLOS Pathogens 17(4): e1009431.](https://doi.org/10.1371/journal.ppat.1009431)


