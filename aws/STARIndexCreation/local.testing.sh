
grabnode 
## 4 cores, gizmof18
## workdir: /fh/scratch/delete90/paguirigan_a/STAR

###gtf=s3://fh-ctr-public-reference-data/genome_data/human/hg38/hg38_refseq_ucsc.gtf
###fasta=s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta
###testData=s3://fh-pi-paguirigan-a/tg/SR/ngs/illumina/apaguiri/150219_SN367_0497_AHGLYGADXX/Unaligned/Project_apaguiri/Sample_NC-23B/
readLength=50

module load awscli
aws s3 cp s3://fh-ctr-public-reference-data/genome_data/human/hg38/Homo_sapiens_assembly38.fasta .
ref_fasta=Gencode_GRCh38.primary_assembly.genome.fa

aws s3 cp s3://fh-ctr-public-reference-data/genome_data/human/hg38/hg38_refseq_ucsc.gtf .
ref_gtf=gencode.v31.annotation.gtf

mkdir rawData
aws s3 sync s3://fh-pi-paguirigan-a/tg/SR/ngs/illumina/apaguiri/150219_SN367_0497_AHGLYGADXX/Unaligned/Project_apaguiri/Sample_NC-23B/ rawData


#Go edit starGenerate.sh to be specific to the species you want to create a reference for.
sbatch starGenerate-human.sh

##Example code for UCSC
cd /fh/scratch/delete90/paguirigan_a/STAR/human/
sbatch --wrap="cd STAR2.7.1_hg38/ && tar -zcf ../STAR2.7.1_hg38.tgz . && cd - " 
sbatch --wrap="cd STAR2.7.1_dm6/ && tar -zcf ../STAR2.7.1_dm6.tgz . && cd - " 
sbatch --wrap="cd STAR2.7.1_canFam3/ && tar -zcf ../STAR2.7.1_canFam3.tgz . && cd - " 
sbatch --wrap="cd STAR2.7.1_mm10/ && tar -zcf ../STAR2.7.1_mm10.tgz . && cd - " 
# This will only tar the files, NOT the whole dir so it will un-tar into a directory you choose in the docker container
## Use Motuz to stick this up into S3 again!



##EXample code for ENSEMBL
sbatch --wrap="cd STAR2.7.1_genomeGRCh38.gencode.v31/ && tar -zcf ../STAR2.7.1_genomeGRCh38.gencode.v31.tgz . && cd - " # This will only tar the files, NOT the whole dir
#TEST THIS sbatch --wrap="tar -xzf STAR2.7.1_genomeGRCh38.tar.gz -C genomeDir"
sbatch --wrap="aws s3 cp STAR2.7.1_genomeGRCh38.gencode.v31.tgz s3://fh-ctr-public-reference-data/genome_data/human/hg38/STAR2.7.1_genomeGRCh38.gencode.v31.tgz"

## DATA!!!  Selected fastq's
r1fastq=/fh/scratch/delete90/paguirigan_a/STAR/rawData/NC-23B_TTAGGC_L001_R1_001.fastq.gz
r2fastq=/fh/scratch/delete90/paguirigan_a/STAR/rawData/NC-23B_TTAGGC_L001_R2_001.fastq.gz

mkdir -p /fh/scratch/delete90/paguirigan_a/STAR/2passalign
cd /fh/scratch/delete90/paguirigan_a/STAR/2passalign

STAR \
    --genomeDir /fh/scratch/delete90/paguirigan_a/STAR/STAR2.7.1_genomeGRCh38.gencode.v31 \
    --readFilesIn ${r1fastq} ${r2fastq} \
    --runThreadN 6 \
    --readFilesCommand zcat \
    --sjdbOverhang 100 \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --quantTranscriptomeBAMcompression 5 

## For RNASeQC
# Running this.
module load Python/3.6.5-foss-2016b-fh3
pip install --user bx-python
python collapse_annotation.py gencode.v31.annotation.gtf gencode.v31.genes.gtf


