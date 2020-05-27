#!/bin/bash
#SBATCH -n 6
#SBATCh --mem 16G


ref_fasta=/fh/scratch/delete90/paguirigan_a/STAR/Gencode_GRCh38.primary_assembly.genome.fa
ref_gtf=/fh/scratch/delete90/paguirigan_a/STAR/gencode.v31.annotation.gtf

module load STAR/2.7.1a-foss-2016b

mkdir -p /fh/scratch/delete90/paguirigan_a/STAR/STAR2.7.1_genomeGRCh38

STAR \
  --runMode genomeGenerate \
  --genomeDir /fh/scratch/delete90/paguirigan_a/STAR/STAR2.7.1_genomeGRCh38 \
  --sjdbGTFfile ${ref_gtf} \
  --genomeFastaFiles ${ref_fasta} \
  --runThreadN 6

