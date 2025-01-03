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
* Levenshtein

## Usage

The pipeline is executed by running the `pacbio-pipeline-with-blast.sh` script. The script requires a FASTA file and a configuration file providing run parameters.

Inputs are provided in the following way:

`./pacbio-pipeline-with-blast.sh $folder $file $config`, where `$folder` is the folder containing the FASTA file, `$file` is the name of the FASTA file (without the file extension), and `$config` is the path to a configuration file.

Upon completion, final single-copy sequences are placed in a folder named `sgs` within the folder containing the specified FASTA input.

## Pipeline Components

`pacbio-pipeline-with-blast.sh` performs a complete analysis towards SGS, although the pipeline is broken down into components and scripts that may be run individually. Major steps are listed below; check the source code for additional information on the syntax of specific scripts.

### Preprocessing: Orientation and Primer Trimming
- `vsearch` is used to orient reads against the provided reference sequence.
- `cutadapt` is used to trim PCR forward and reverse primers. 
	
### UMI Binning
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
- UMIs that are 1-base edit distance away from each other are compared; if one bin is greater than 2x the size of the other, the smaller bin is deemed fake. Its read count is added to the larger bin, but the bin is not re-processed.
- PCR errors for each bin are compiled; error positions are used as fingerprints to associate bins with each other. UMIs with 1-base edit distance that share fingerprints are merged into the largest bin regardless of read count imbalance.
- A mixture model of real vs. fake bins is optimized by the bin size data (real: Gaussian; fake: exponential) and any known fake UMIs. Any bins with Prob(fake | bin size) ≥ 0.5 are discarded.

### Final SGS Alignment
- Polished consensus sequences are written to the `$folder/sgs/` folder. Sequences are provided in 3 schemes: intra-sample aligned, aligned against the reference, and unaligned.

## Post Processing of SGS

Before proceeding, it is highly recommended to manually review and curate the alignment. Truncated and/or defective sequences may be removed at this stage if they are encountered. The alignment of all insertions, deletions, and homopolymer regions should be reviewed.

### Haplotype Calling
- To call haplotypes, run `revert-mutations-reference.py` on a curated alignment of the called single-genome sequences.
- Example: `python3 scripts/revert-mutations-reference.py path/to/sgs-alignment.fasta`
- This step is not automatically performed after calling single-genome sequences because it is highly recommended to manually review the final SGS alignment before post-processing sequences.
- A high-quality alignment is critical to ensure proper variant calling and minimize false positives.
- Haplotypes will only be called within the region that is included as a reference (the reference is assumed to be the first sequence in the curated alignment).

## Citation

If you use this pipeline in your work, please cite:

[Ko SH, Radecki P, et al. (2024) Rapid intra-host diversification and evolution of SARS-CoV-2 in advanced HIV infection. Nature Communications 15(7240). doi: 10.1038/s41467-024-51539-8](https://doi.org/10.1038/s41467-024-51539-8)

