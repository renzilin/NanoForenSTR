# NanoForenSTR v0.0.1

This is the script for Forensic STRs genotyping in the manuscript: *Nanopore sequencing of forensic STRs and SNPs using Verogen’s ForenSeq DNA Signature Prep Kit and MinION.* We've tested it on 54 STRs (27 autosomal STRs, 7 X-STRs, and 20 Y-STRs). In details, 22 of 54 STRs can be genotyped robustly and correctly in 30 collected samples and 3 standard 2800 M duplicates.



## Workflow

In our manuscript's workflow, (1) nanopore reads were aligned to the human reference genome (hg19) using Minimap2 and transformed to the `*.bam` file using Samtools; (2) using the `NanoForenSTR.py`, with the required `*.bam` file and pattern file (contains 5 columns (STR name, chromosome, start, end, repeat unit)), to genotype STRs. In the folder `test/pattern_file`, we provide a pattern file including 54 STRs used in our manuscript.

NanoForenSTR firstly extracts reads (the default flanking length is 30 bp) from *.bam* file according to the STR coordinate from pattern file in the genotyping process. Secondly, the repeat units in the pattern file are searched in each extracted read. Then, the read is divided into several fragments, some of which are composed to repeat units. For the rest fragments that are not composed to repeat units, we use the Smith-Waterman (SW) algorithm to align the fragments against a sequence made up of repeat units and has the same length as the fragment. If the alignment result shows less than one mismatch or two gaps in each repeat unit part, this fragment will be determined as a part of the repeat region. After SW alignment process, the repeat regions are determined, and the copy number of repeat units in this region can be counted. Finally, we make a summary of the copy number from reads. The first allele is the copy number with the largest number of supporting reads. The second allele is the copy number whose number of supporting reads is over 50% of the first allele. Otherwise, this is a homozygous allele. 



## Installation

We recommend to use  `conda` to build NanoForenSTR environment and run our script under it. We've tested it in Ubuntu platform and Windows Subsystem Linux 2 (WSL2) platform.

1.  Build the environment:

```bash
## create a conda environment
conda create -n NFS-v0.0.1
conda activate NFS-v0.0.1

## install dependent packages
conda install -c conda-forge -c bioconda  python=3.8 pysam=0.16 numpy=1.19 pandas=1.1.3 cython=0.29 tqdm

cd myLibs

python setup_myDPAlign.py build_ext --inplace

```



2. Download the code from our repo and test

```bash
git clone https://github.com/renzilin/NanoForenSTR.git
cd NanoForenSTR/test

python ../NanoForenRepeat.py LA --BAM bam_file/barcode01.test.bam --PAT pattern_file/STRtest.pat --ID test
```



## Usage

In our manuscript, we use the option `LA`, which is also recommended, to genotype the Forensic STRs. There are 2 options for genotyping STRs. The first one is quick mode, which is a test version. The second mode uses SW to handle the mismatch and gap. 

### Quick Start

 ```bash
conda activate NFS-v0.0.1

NanoForenSTR="[PATH_TO_NanoForenSTR]"

## required arguments
BAM_FILE="[PATH_TO_BAM_FILE]"
PATTERN_FILE="[PATH_TO_PAT_FILE]"
SAMPLE_NAME="[SAMPLE_NAME]"

python $NanoForenSTR/NanoForenSTR.py LA \
--BAM $BAM_FILE \
--PAT $NanoForenSTR/test/pattern_file/STR.pat \
--ID  $SAMPLE_NAME
 ```



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







