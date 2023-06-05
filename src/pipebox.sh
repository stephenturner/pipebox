#!/bin/bash
set -e

# Check usage
if [ "$#" -ne 4 ]; then echo "Usage: pipebox.sh <read1> <read2> <reffa> <outbase>"; exit 1; fi

# Set variables used throughout this script
srcdir=$(dirname $0)
read1=$1
read2=$2
reffa=$3
outbase=$4

# Check that files exist inside the container. 
if [ ! -f "$read1" ]; then echo "$read1 doesn't exist"; exit 1; fi
if [ ! -f "$read2" ]; then echo "$read2 doesn't exist"; exit 1; fi
if [ ! -f "$reffa" ]; then echo "$reffa doesn't exist"; exit 1; fi

# Check that a samtools faidx index exists. If not, create it.
if [ ! -f "$reffa.fai" ]; then samtools faidx $reffa; fi

# Check that a bwa index exists. If not, create it.
if [ ! -f "$reffa.bwt" ]; then bwa index $reffa; fi

# Make a directory for the results if it doesn't yet exist
mkdir -p $(dirname $outbase)

# Read QC
echo "Read QC"
seqtk mergepe $read1 $read2 | seqtk fqchk - > $outbase.seqtkfqchk.tsv

# Alignment
echo "Alignment"
rg="@RG\tID:${outbase}\tSM:${outbase}\tLB:${outbase}"
bwa mem -R $rg $reffa $read1 $read2 | samtools sort -Obam -o $outbase.bam
samtools stats $outbase.bam > $outbase.samtoolsstats.tsv

# Variant calling
echo "Variant calling"
bcftools mpileup -f $reffa $outbase.bam | bcftools call -mv - | bcftools sort -Oz -o $outbase.vcf.gz
bcftools stats $outbase.vcf.gz > $outbase.bcfstats.tsv

# Run the R script for post-processing results
echo "Postprocessing results"
Rscript $srcdir/pipebox-post.R $outbase

echo "Done! See results at $outbase*"