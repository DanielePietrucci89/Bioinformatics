#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=20 --mem=110GB
#SBATCH --time=18:00:00
#SBATCH --job-name QC
#SBATCH --account=
#SBATCH --partition=g100_usr_prod
#SBATCH --mail-user=
#SBATCH	--mail-type=ALL
#SBATCH --output=QC.out
#SBATCH --error=QC.err

###################################################################################
# This script is used to perform quality control analyses on			  #
# aDNA sequencing data related to medieval ceramics                               #
# All data must be placed in a folder named 000_RawData                           #
###################################################################################

################# PATH & VARIABLES ################################################
RAW_DATA="000_RawData"
FASTQC_1="001_Fastqc/000_RawData"
FASTQC_2="001_Fastqc/001_PostCutadapt"
FASTQC_3="001_Fastqc/002_PostTrimming"
CUTADAPT="002_Quality_Filtered/000_Cutadapt"
TRIMMOMATIC="002_Quality_Filtered/001_Trimmomatic"
MERGED="003_Merged"
###################################################################################

###################################################################################
module load profile/advanced
module load profile/bioinf
module load autoload fastqc/0.11.9
module load autoload trimmomatic/0.39

### ---> FIRST STEP: Evaluating raw data quality using fastqc
mkdir -p $FASTQC_1
#fastqc -t 4 -o $FASTQC_1 -f fastq $RAW_DATA/*
# Inspect the results

### ---> SECOND STEP: Remove the adapter sequence using Cutadapt
mkdir -p $CUTADAPT
source activate /../../../../miniconda3/envs/cutadaptenv
for forward in `ls $RAW_DATA/*_R1_001.fastq.gz`; do
  reverse=${forward/_R1_/_R2_}
  ls -lhrt $forward
  ls -lhrt $reverse
	name_for=$(basename $forward)
  name_rev=$(basename $reverse)
	out_for=$CUTADAPT/$name_for
	out_rev=$CUTADAPT/$name_rev
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 16 -o $out_for -p $out_rev $forward $reverse
done
conda deactivate

mkdir -p $FASTQC_2
#fastqc -t 4 -o $FASTQC_2 -f fastq $CUTADAPT/*

### ---> THIRD STEP: Remove low-quality reads
mkdir -p $TRIMMOMATIC
for forward in `ls $CUTADAPT/*R1*.fastq.gz`; do
	reverse=${forward/R1/R2}
	ls $forward
	ls $reverse
	basename_file=$(basename $forward)
	output_trimlog=$TRIMMOMATIC/${basename_file/_L001_R1_001.fastq.gz/_trimlog.txt}
	output_for_pai=$TRIMMOMATIC/${basename_file/_L001_R1_001.fastq.gz/_R1_pai.fastq.gz}
	output_rev_pai=$TRIMMOMATIC/${basename_file/_L001_R1_001.fastq.gz/_R2_pai.fastq.gz}
	output_for_unp=$TRIMMOMATIC/${basename_file/_L001_R1_001.fastq.gz/_R1_unp.fastq.gz}
	output_rev_unp=$TRIMMOMATIC/${basename_file/_L001_R1_001.fastq.gz/_R2_unp.fastq.gz}
	java -jar $TRIMMOMATIC_HOME/bin/trimmomatic-0.39.jar PE -phred33 -threads 12 \
		-trimlog $output_trimlog \
		$forward $reverse \
		$output_for_pai $output_for_unp \
		$output_rev_pai $output_rev_unp \
		MINLEN:30 AVGQUAL:30
done

mkdir -p $FASTQC_3
#fastqc -t 4 -o $FASTQC_3 -f fastq $TRIMMOMATIC/*

### --> FOURTH STEP: Merging using LeeHom
mkdir -p $MERGED
source activate /../../../../miniconda3/envs/leehomEnv
for forward in `ls $TRIMMOMATIC/*R1_pai*.fastq.gz`; do
        reverse=${forward/R1/R2}
        ls $forward
        ls $reverse
        basename_file=$(basename $forward)
        basename_file=${basename_file/_R1_pai.fastq.gz/.fastq.gz}
	out_filename=$MERGED/$basename_file
        echo $out_filename
	leeHom -fq1 $forward -fq2 $reverse -fqo $out_filename --ancientdna 
done
conda deactivate
