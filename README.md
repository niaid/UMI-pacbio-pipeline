# PacBio UMI Analysis Pipeline

A command-line pipeline for analysis of PacBio sequencing data from libraries with unique molecular identifiers (UMIs).

Inputs to the pipeline are FASTA PacBio CCS reads and a configuration file describing necessary parameters, such as a reference file and primer sequences.

Outputs include UMI-consensus target sequences (i.e., single-copy/single-genome/single-molecule sequences) as well as intermediate outputs from various processing steps.  

The pipeline assumes a single target sequence and the following read structure:

*5'* – **[ PCR forward primer ][ target sequence ][ RT primer ][ UMI ][ PCR reverse primer ]** – *3'*

The primer sequences and UMI region can be of any length. Reads need not be oriented prior to running the pipeline; the pipeline will orient the reads against the provided target sequence reference.

*This pipeline was developed by the Virus Persistence and Dynamics Section at the National Institutes of Health (NIH).*


## Requirements

This pipeline requires several command-line tools as well as a Python 3 environment.

The following command-line tools are required:
* [bcftools and samtools](https://www.htslib.org/download/)
* [BLAST+ / blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [mafft](https://mafft.cbrc.jp/alignment/software/)
* [minimap2](https://github.com/lh3/minimap2)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
* [vsearch](https://github.com/torognes/vsearch)

The following Python packages are also required:
* Biopython
* joblib
* tqdm
* networkx
* numpy
* pandas
* pysam
* scipy
* matplotlib

## Usage

The pipeline is executed by running the `pacbio-pipeline-with-blast.sh` script. The script requires a FASTA file and a configuration file providing run parameters.

Inputs are provided in the following way:

`./pacbio-pipeline-with-blast.sh $folder $file $config`, where `$folder` is the folder containing the FASTA file, `$file` is the name of the FASTA file (without the file extension), and `$config` is the path to a configuration file.

Upon completion, final single-copy sequences are placed in a folder named `sgs` within the folder containing the specified FASTA input.

## Pipeline components:

`pacbio-pipeline-with-blast.sh` performs a complete analysis, although the pipeline is broken down into components and scripts that may be run individually. Major steps are listed below; check the source code for additional information on the syntax of specific scripts.

### Preprocessing: Orientation and Primer Trimming
- `vsearch` is used to orient reads against the provided reference sequence.
- `cutadapt` is used to trim PCR forward and reverse primers. 
	
### UMI Binning (`bin_umis.py`)
- Attempts to match RT primer sequence to trimmed reads; if successful, the UMI is extracted and stored.
- The relevant section of the read can differ from the provided RT primer sequence by up to 1 base.
- UMIs with 10 CCS or more are kept as putative bins for consensus calling.

### UMI Consensus Determination
- `vsearch` is used to cluster reads within bins; only bins with 1 coherent cluster (largest cluster is at least 2x larger than next cluster *and* largest cluster has at least 10 reads) are kept.
- Reads are mapped with `minimap2` to the bin consensus determined by `vsearch`.
- `bcftools` is used to call variants from the read mapping and determine the final consensus sequence.

### Sequence Filtering
- Bin consensus sequences are compared to the reference sequence with `blastn`; only sequences which are similar to the reference are kept.

### Fake UMI Removal
- UMIs that are 1-base edit distance away from each other are compared; if one bin is greater than 2x the size of the other, the smaller bin is deemed fake. Its read count is added to the larger bin, but the reads are not re-processed.
- PCR errors for each bin are compiled; error positions are used as fingerprints to associate bins with each other. UMIs with 1-base edit distance that share fingerprints are merged.

### PCR Error Correction
- PCR error positions (ambiguous bases in a bin) are checked and reverted to the consensus if the position is highly conserved in the sample (>90% identity) and one of the conflicting nucleotides matches the sample consensus. 

### Final Alignment
- Polished consensus sequences are written to the `$folder/sgs/` folder. Sequences are provided in 3 schemes: intra-sample aligned, aligned against the reference, and unaligned.


## Citation

If you use this pipeline in your work, please cite:

[Ko SH, Bayat Mokhtari E, Mudvari P, Stein S, Stringham CD, et al. (2021) High-throughput, single-copy sequencing reveals SARS-CoV-2 spike variants coincident with mounting humoral immunity during acute COVID-19. PLOS Pathogens 17(4): e1009431.](https://doi.org/10.1371/journal.ppat.1009431)


