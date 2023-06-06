#!/bin/bash
set -e
srcdir=$(dirname $0)
read1=$1
read2=$2
reffa=$3
outbase=$4
if [ ! -f "$reffa.fai" ]; then samtools faidx $reffa; fi
if [ ! -f "$reffa.bwt" ]; then bwa index $reffa; fi
mkdir -p $(dirname $outbase)
seqtk mergepe $read1 $read2 | seqtk fqchk - > $outbase.seqtkfqchk.tsv
rg="@RG\tID:${outbase}\tSM:${outbase}\tLB:${outbase}"
bwa mem -R $rg $reffa $read1 $read2 | samtools sort -Obam -o $outbase.bam
samtools stats $outbase.bam > $outbase.samtoolsstats.tsv
bcftools mpileup -f $reffa $outbase.bam | bcftools call -mv - | bcftools sort -Oz -o $outbase.vcf.gz
bcftools stats $outbase.vcf.gz > $outbase.bcfstats.tsv
Rscript $srcdir/pipebox-post.R $outbase