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

module load profile/advanced
module load profile/bioinf
module load autoload bwa/0.7.12
module load autoload samtools/1.13
#module load autoload kraken/2

################# PATH & VARIABILES ################################################

# Read here for the DB instructions!
# Download standard db (from https://benlangmead.github.io/aws-indexes/k2)
# Downloading december 2024 STANDARD database because 2025 DB has issue with lineage notations (domain is R2 in 2025) 
# wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20241228.tar.gz (it's a huge database, it may requires a LOT of time!)
# once dowloaded decompress the file with tar -xvzf,
# then, the Kraken2 DB and Bracken2 DB were stored inside the path of $DB variable

ALN_OUT_DIR_Hum_1="005_ALN_OUT_DIR_HUMAN/000_ALN"
RAW_DATA_BACTERIA="008_Raw_Data_Only_Bacteria"
ALN_OUT_DIR_K2_MP="009_ALN_K2_MP"
ALN_OUT_DIR_KRAKEN="009_ALN_K2_MP/000_KRAKEN_RESULTS"
ALN_OUT_DIR_METAPHLAN="009_ALN_K2_MP/001_METAPHLAN_RESULTS"
KRAKEN_DB="/../../../..l/kraken_db"

####################################################################################

### Extract reads which were not aligned on the human genome 

mkdir -p $RAW_DATA_BACTERIA
for input_file in `ls $ALN_OUT_DIR_Hum_1/*.sorted.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/.sorted.bam/}
        NAME=$file_basename
        samtools fastq -f 4 $input_file | gzip > $RAW_DATA_BACTERIA/$output_basename".fastq.gz"
done

### Run Kraken2
### Reads that were not aligned on the human genome were aligned using Kraken2
### This quick analysis allows the investigation of potential microbial taxa which may be present in the sample
### Interesting taxa (i.e. bacteria that may be involved in human pathology or commensal bacteria)
###   will be further investigated in another script, after checking the results of Kraken2

source activate /../../../../miniconda3/envs/Kraken-Bracken
mkdir -p $ALN_OUT_DIR_K2_MP
mkdir -p $ALN_OUT_DIR_KRAKEN
for input_file in `ls $RAW_DATA_BACTERIA/*.fastq.gz`; do
	file_basename=$(basename $input_file)
	output_basename=${file_basename/.fastq.gz/}
	KRAKEN_REPORT=$ALN_OUT_DIR_KRAKEN/$output_basename".report"
	KRAKEN_OUTPUT=$ALN_OUT_DIR_KRAKEN/$output_basename".kraken"
  SPECIES_MPA=$ALN_OUT_DIR_KRAKEN/$output_basename".species_only.mpa.report"
	BACTERIA_MPA=$ALN_OUT_DIR_KRAKEN/$output_basename".species_bacteria_only.mpa.report"

	# All taxa
	echo "starting Kraken2 for $file_basename"
	kraken2 --db $KRAKEN_DB --threads 18 --report $KRAKEN_REPORT --output $KRAKEN_OUTPUT --use-mpa-style --gzip-compressed $input_file

	# Only bacteria
	# Filter species-only assignments from MPA-style report
  grep "s__" $KRAKEN_REPORT > $SPECIES_MPA
	grep -E "d__Archaea|d__Bacteria" $SPECIES_MPA > $BACTERIA_MPA
done
conda deactivate

#exit


#merging all three reports into one file
#filtered bacteria reports
python /g100_work/ELIX6_gruber/TOOLS/KrakenTools/combine_mpa.py \
  -i $ALN_OUT_DIR_KRAKEN/48064_ID1564_12_1A_S12.species_bacteria_only.mpa.report $ALN_OUT_DIR_KRAKEN/48065_ID1564_13_1C_S13.species_bacteria_only.mpa.report \
     $ALN_OUT_DIR_KRAKEN/48066_ID1564_14_14A_S14.species_bacteria_only.mpa.report $ALN_OUT_DIR_KRAKEN/48067_ID1564_15_14C_S15.species_bacteria_only.mpa.report \
     $ALN_OUT_DIR_KRAKEN/48068_ID1564_16_20A_S16.species_bacteria_only.mpa.report $ALN_OUT_DIR_KRAKEN/48069_ID1564_17_20B_S17.species_bacteria_only.mpa.report \
     $ALN_OUT_DIR_KRAKEN/48070_ID1564_18_20C_S18.species_bacteria_only.mpa.report \
  -o $ALN_OUT_DIR_K2_MP/kraken2_merged_reports_species_bacteria_only.txt

# Modify column names, following the name of samples provided in input
sed -i $'1s/.*/#Classification\t#1A\t#1C\t#14A\t#14B\t#20A\t#20B\t#20C/' $ALN_OUT_DIR_K2_MP/kraken2_merged_reports_species_bacteria_only.txt

exit
