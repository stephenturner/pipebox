#!/bin bash

# Illumina sequencing data from a yeast genome sequencing project from CSHL
# https://www.ncbi.nlm.nih.gov/bioproject/167999

# Get 40,000 reads
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507778/SRR507778_1.fastq.gz | gzip -d | head -n 40000 | gzip -c > SRR507778_1.fastq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR507/SRR507778/SRR507778_2.fastq.gz | gzip -d | head -n 40000 | gzip -c > SRR507778_2.fastq.gz

# get the yeast reference genome
curl ftp://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz > yeastref.fa.gz

# Decompress it
gunzip -k yeastref.fa.gz

# Index it
samtools faidx yeastref.fa
bwa index yeastref.fa

