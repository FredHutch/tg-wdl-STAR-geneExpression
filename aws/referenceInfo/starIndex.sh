#!/bin/bash
#SBATCH -n 6 		# number of nodes requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem 16G
#SBATCH --job-name STAR_index 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written

mkdir -p /fh/scratch/delete90/paguirigan_a/STAR
WORKDIR=/fh/scratch/delete90/paguirigan_a/STAR

module load STAR/2.7.1a-foss-2016b

mkdir -p ${WORKDIR}/STAR_2.7.1a_GRCh28_gencode_index_50base

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${WORKDIR}/STAR_2.7.1a_GRCh28_gencode_index_50base \
--genomeFastaFiles /fh/fast/paguirigan_a/trgen/ReferenceData/hg38/Gencode_GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /fh/fast/paguirigan_a/trgen/ReferenceData/hg38/gencode.v29.annotation.gtf \
--sjdbOverhang 49 # For 50 base PE sequencing