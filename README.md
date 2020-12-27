# NanoForenSTR

The script for genotyping Forensic STR loci in *Nanopore sequencing of forensic STRs and SNPs using Verogen’s ForenSeq DNA Signature Prep Kit and MinION.*



## Features

1. We genotyped XX autosomal STRs, xx X-STRs, and xx Y-STRs using our script. 
2. The whole working process can be explained by the figures below.





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

In our 







