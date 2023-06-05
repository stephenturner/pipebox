seqtk mergepe mysamp_R1.fastq.gz  mysamp_R2.fastq.gz | seqtk fqchk - > mysamp.seqtkfqchk.tsv
bwa mem -R '@RG\tID:PRJNA167999\tSM:SRR507778\tLB:SRX151723' myref.fasta mysamp_R1.fastq.gz mysamp_R2.fastq.gz | samtools sort -O bam -o mysamp.bam
samtools stats mysamp.bam > mysamp.samtoolsstats.tsv
bcftools mpileup -f myref.fasta mysamp.bam | bcftools call -mv - | bcftools sort -Oz -o mysamp.vcf.gz
bcftools stats mysamp.vcf.gz > mysamp.bcfstats.tsv
