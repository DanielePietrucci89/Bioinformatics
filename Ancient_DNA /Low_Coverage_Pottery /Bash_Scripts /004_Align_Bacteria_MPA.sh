#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=32 --mem=150GB
#SBATCH --time=10:00:00
#SBATCH --job-name aln_H
#SBATCH --account=
#SBATCH --partition=g100_usr_prod
#SBATCH --mail-user=
#SBATCH	--mail-type=ALL
#SBATCH --output=aln_H.out
#SBATCH --error=aln_H.err

####################################################################################

####################################################################################

# NOTES

# Install metaphlan if you never did!
# conda create -n metaphlan -c bioconda -c conda-forge metaphlan bowtie2 -y
# this is metaphlan version 4.2.4 (21 Oct 2025) 
# dowload the database using the command:
# metaphlan --install --db_dir /g100_work/IscrB_VIROME/MetaPhlAn/MetaPhlAn_DB
# DB should be place in a directory defined by the user, the same directory should be in METAPHLAN_DB Variable

####################################################################################

module load profile/advanced
module load profile/bioinf

################# PATH & VARIABILES ################################################

ALN_OUT_DIR_Hum_1="005_ALN_OUT_DIR_HUMAN/000_ALN"
RAW_DATA_BACTERIA="008_Raw_Data_Only_Bacteria"
ALN_OUT_DIR_K2_MP="009_ALN_K2_MP"
ALN_OUT_DIR_METAPHLAN="009_ALN_K2_MP/001_METAPHLAN_RESULTS"
METAPHLAN_DB="/../../MetaPhlAn/MetaPhlAn_DB"

####################################################################################

### Extract reads which were not aligned on the human genome 
### If you already did it for Kraken 2 skip these lines
mkdir -p $RAW_DATA_BACTERIA
for input_file in `ls $ALN_OUT_DIR_Hum_1/*.sorted.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/.sorted.bam/}
        NAME=$file_basename
        samtools fastq -f 4 $input_file | gzip > $RAW_DATA_BACTERIA/$output_basename".fastq.gz"
done

source activate /../../../../miniconda3/envs/metaphlan
mkdir -p $ALN_OUT_DIR_K2_MP
mkdir -p $ALN_OUT_DIR_METAPHLAN
for input_file in `ls $RAW_DATA_BACTERIA/*.fastq.gz`; do
	file_basename=$(basename $input_file)
	output_basename=${file_basename/.fastq.gz/}
	MPA_MAPOUT=$ALN_OUT_DIR_METAPHLAN/$output_basename".bowtie2.bz2"
	MPA_OUTPUT=$ALN_OUT_DIR_METAPHLAN/$output_basename"_profile.mpa"
  echo "Processing $input_file..."
  metaphlan $input_file \
  --input_type fastq \
  --db_dir $METAPHLAN_DB \
  --mapout $MPA_MAPOUT \
  --nproc 16 \
  -o $MPA_OUTPUT
  echo "Finished $input_file"
done

merge_metaphlan_tables.py $ALN_OUT_DIR_METAPHLAN/*"_profile.mpa" > $ALN_OUT_DIR_K2_MP/merged_metaphlan_report.tsv

conda deactivate

exit
