#!/bin/bash
#SBATCH -n 6
#SBATCh --mem 16G

module purge

species=drosophila
version=dm6
mkdir -p /fh/scratch/delete90/paguirigan_a/STAR/${species}
mkdir -p /fh/scratch/delete90/paguirigan_a/STAR/${species}/STAR2.7.1_${version}

genome=/fh/scratch/delete90/paguirigan_a/STAR/${species}/genome.fa
genes=/fh/scratch/delete90/paguirigan_a/STAR/${species}/genes.gtf

module load awscli

aws s3 cp s3://fh-ctr-public-reference-data/genome_data/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa ${genome}
aws s3 cp s3://fh-ctr-public-reference-data/genome_data/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa.fai ${genome}.fai
aws s3 cp s3://fh-ctr-public-reference-data/genome_data/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf ${genes}

module load STAR/2.7.1a-foss-2016b

STAR \
  --runMode genomeGenerate \
  --genomeDir /fh/scratch/delete90/paguirigan_a/STAR/${species}/STAR2.7.1_${version} \
  --sjdbGTFfile ${genes} \
  --genomeFastaFiles ${genome} \
  --runThreadN 6

## 100bp read length is default
