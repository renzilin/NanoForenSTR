# NanoForenSTR

This script is for Forensic STRs genotyping in the manuscript: *Nanopore sequencing of forensic STRs and SNPs using Verogen’s ForenSeq DNA Signature Prep Kit and MinION.* We've tested it on 54 STRs (XX autosomal STRs, xx X-STRs, and xx Y-STRs). In details, XX of XX STRs can be genotyped robustly and correctly in 30 collected samples and 3 standard 2800 M duplicates.



## Workflow

In our manuscript's workflow, (1) nanopore reads were aligned to the human reference genome (hg19) using Minimap2 and transformed to the `*.bam` file using Samtools; (2) using the `x.py` to genotype STRs with the required `*.bam` file and pattern file that contains 5 columns (STR name, chromosome, start, end, repeat unit). In the folder `x/x`, we provide a pattern file including 54 STRs used in our manuscript.

In genotype process, NanoForenSTR firstly extract reads from *.bam* file according to the STR coordinate from pattern file. Secondly, the repeat unit in pattern file are searched in each extracted read. Then, the read is divided into several fragments, some of which are composed of repeat units, and some of which are not. As for the fragments which are not composed of repeat units, we use the Smith-Waterman algorithm to align the fragments against a sequence of the same length as the fragment, which is made up of repeat units. If the alignment result shows less than 1 mismatch or 2 gaps in each repeat unit part, this fragment will be determined as a part of repeat region. After SW alignment process, the repeat region are determined and the copy number of repeat unit in this region can be counted. Finally, we make a summary of copy number from reads. The first allele is the copy number with the largest number of supporting reads. The second allele is the copy number whose number of supporting reads is at least 50% of that of the first allele. Otherwise, this is a homozygous allele. 



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

In our manuscript, we use the option `LA`, which is also recommended, to genotype the Forensic STRs. There are 2 options for genotyping STRs. The first one is quick mode, which is just a test version. The second mode uses SW to handle the mismatch and gap problem 





### General Usage Options

```bash
usage: python NanoForenRepeat.py <options> [<args>]

Available options are:    
    quick    Genotyping forensic STR by quick mode
    LA       Genotyping forensic STR by local align (recommend)
    
    
    
Genotyping forensic STR

positional arguments:
  {quick,LA}  sub-command help
    quick     Genotyping forensic STR by quick mode
    LA        Genotyping forensic STR by local align

optional arguments:
  -h, --help  show this help message and exit

```







