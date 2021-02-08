# NanoForenSTR

This script is for Forensic STRs genotyping in the manuscript: *Nanopore sequencing of forensic STRs and SNPs using Verogen’s ForenSeq DNA Signature Prep Kit and MinION.* We've tested it on 54 STRs (XX autosomal STRs, xx X-STRs, and xx Y-STRs). In details, XX of XX STRs can be genotyped robustly and correctly in 30 collected samples and 2800 M. ()



## Workflow

In our manuscript's workflow, (1) long reads were aligned to the human reference genome (hg19) using Minimap2 and transform to the `*.bam` file using samtools; (2) using the `x.py` to genotype STRs with the required `*.bam` file and pattern file that contains 5 columns (STR name, chromosome, start, end, repeat unit). In the folder `x/x`, we provide a pattern file that used in our manuscript including 54 STRs.

In genotype process, NanoForenSTR firstly search the repeat unit provided by pattern file in aligned reads. Then, the read is split into several sub-regions which is . For the mismatching sub-regions, NanoForenSTR generates  

we use the Smith-Waterman algorithm to align the sub-region against repeat unit tolerating 1 mismatch or 2 gaps ( This is from our testing ). For the mismatch sub-region, if 

If the mismatching sub-region is aligned, then the 



## Installation

We recommend you to use  `conda` to build NanoForenSTR environment and run our script under it. 

1.  Build the environment:

```bash
## create a conda environment
conda create -n NFS-v0.1

## enter the environment
conda activate NFS-v0.1

## install dependent packages
conda install -c conda-forge -c bioconda  python=3.8 pysam=0.16 numpy=1.19 pandas=1.1.3 cython=0.29
```



2. Download the code from our repo:

```bash
git clone https://github.com/renzilin/NanoForenSTR.git
```



3. Convert the *.pyx to *.c for accelaration purpose.

```bash
cd myLibs
python setup_myPairwiseAlignment.py  build_ext --inplace
```







## Usage

In our manuscript, we use the option `LA` to  genotype the Forensic STRs. The 





### General Usage Options

```bash
usage: python ForenRepeat.py <command> [<args>]

Available options are:    
    quick    Repeat quantification by quick mode
    LA       Repeat quantification by local align LA
       [-h] [--BAM BAM] [--PAT PAT] [--ID ID]

Required arguments:
  -h, --help  show this help message and exit
  --BAM BAM   input1: *.bam file
  --PAT PAT   input2: repeat pattern file
  --ID ID     output file name

```







