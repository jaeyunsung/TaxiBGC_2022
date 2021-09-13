# TaxiBGC: Taxonomy-guided Identification of Biosynthetic Gene Clusters

## Description 
TaxiBGC (Taxonomy-guided Identification of Biosynthetic Gene Clusters): A computational pipeline for identifying experimentally verified Biosynthetic Gene Clusters (BGCs) and inferring their annotated secondary metabolites (SMs) from metagenomic shotgun sequencing data. The TaxiBGC pipeline includes three major steps: 1) species-level, taxonomic profiling on the metagenome; 2) a first-pass prediction of BGCs through querying the species (identified in the first step) in the TaxiBGC database; and 3) confirmation (in silico) of the predicted BGCs (from the second step) based on the detection of BGC genes in the metagenome.
MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea, and Eukaryotes) from metagenomic shotgun sequencing data (not 16S) at species-level resolution.

If you use TaxiBGC, please cite:

Bakshi and Gupta et al. TaxiBGC: a Taxonomy-guided Approach for the Identification of Experimentally Verified Microbial Biosynthetic Gene Clusters in Shotgun Metagenomic Data 
DOI to bioRxiv preprint: https://doi.org/10.1101/2021.07.30.454505

-------------

## Pre-requisites
1. R version 3.6.3 (or higher)
2. R library 'dplyr' and 'argparse'
3. MetaPhlAn (https://github.com/biobakery/MetaPhlAn)
4. Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
5. SAMtools (http://www.htslib.org/download/)

## Installation
TaxiBGC can be downloaded directly from GitHub either as a compressed archive (https://github.com/jaeyunsung/TaxiBGC_2021/archive/refs/heads/main.zip) or using the git command line client:

$ git clone https://github.com/jaeyunsung/TaxiBGC_2021.git

This will create a directory 'TaxiBGC_2021-main' in the current working directory, and will include the following:
1. TaxiBGC_script.R: This is the main R script of the TaxiBGC pipeline.
2. database.zip: After unzipping 'database.zip', it will create a directory 'database' inside the 'TaxiBGC_2021-main' directory. This contains the following:
	1. TaxiBGC_database.csv: This is the essential background database for the TaxiBGC pipeline, which contains information about BGCs and their corresponding SMs
	2. A directory 'bgc_fasta', which includes nucleotide fasta sequences for each BGC included in the TaxiBGC background database ('TaxiBGC_database.csv')
3. example: This directory contains the input paired-end metagenome ('example_1.fq.gz' and 'example_2.fq.gz') and output file ('BGCs.csv') for a demo run of TaxiBGC.
4. ReadME.md

First, unzip 'database.zip' found inside the 'TaxiBGC_2021-main' directory. Then, install all the TaxiBGC dependencies, including MetaPhlan, Bowtie2, and SAMtools. The best way to install MetaPhlan, Bowtie2, and SAMtools is through Conda via the Bioconda channel (https://bioconda.github.io/user/install.html). Alternatively, these tools can be manually installed from their source and added to system $PATH.

-------------

## Input files
TaxiBGC accepts paired-end (PE) metagenomic sequences in .fastq or .fastq.gz formats.

## Output file
TaxiBGC generates only one "BGCs.csv" file inside the output directory. This file includes the predicted BGCs and their corresponding SMs for the corresponding metagenome.

## How to run TaxiBGC
You can run TaxiBGC using the following command in the terminal (Linux/Mac):

$ Rscript --vanilla $DIR/TaxiBGC_2021-main/TaxiBGC_script.R -R1 $DIR/FORWARD_END_METAGENOME_READS -R2 $DIR/REVERSE_END_METAGENOME_READS -o $DIR_TO_OUTPUT_FILES -d $DIR/TaxiBGC_2021-main/database

## Demo run of TaxiBGC
$ cd $DIR/TaxiBGC_2021-main

$ Rscript --vanilla TaxiBGC_script.R -h

$ Rscript --vanilla $DIR/TaxiBGC_2021-main/TaxiBGC_script.R -R1 $DIR/TaxiBGC_2021-main/example/example_1.fq.gz -R2 $DIR/TaxiBGC_2021-main/example/example_2.fq.gz -o $DIR/TaxiBGC_2021-main/demo -d $DIR/TaxiBGC_2021-main/database

After completion of the TaxiBGC run on the example metagenome (example_1.fq & example_2.fq), an output file 'BGCs.csv' will be generated inside $DIR/TaxiBGC_2021-main/demo.
