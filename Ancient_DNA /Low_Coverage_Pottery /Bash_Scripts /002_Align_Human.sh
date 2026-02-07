#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=20 --mem=110GB
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

####################################################################################

module load profile/advanced
module load profile/bioinf
module load autoload bwa/0.7.12
module load autoload samtools/1.13
module load autoload picard/2.25.7

################# PATH & VARIABILES ################################################

MERGED="003_Merged/"
DB="004_DB_HUMAN"
REF=$DB/"genome.fa"
ALN_OUT_DIR_1="005_ALN_OUT_DIR_HUMAN/000_ALN"
ALN_OUT_DIR_2="005_ALN_OUT_DIR_HUMAN/001_PICARD"
ALN_OUT_DIR_3="005_ALN_OUT_DIR_HUMAN/002_QUAL"
ALN_OUT_DIR_4="005_ALN_OUT_DIR_HUMAN/003_PRIMARY"

NCORES=16
SEED=1000
GAP=2
EDIT=0.1
DICT="004_DB/all_RefSeq.dict"


####################################################################################

# Download the Human Genome from ENSEMBL, decompress the genome and use bwa index
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mkdir -p $DB
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz $DB/genome.fa.gz
gunzip -d $DB/genome.fa.gz
bwa index $DB/genome.fa

echo "DB used: $REF"

####################################################################################

mkdir -p $ALN_OUT_DIR_1
for input_file in `ls $MERGED/*.fastq.gz.fq.gz`; do
	file_basename=$(basename $input_file)
	output_basename=${file_basename/.fastq.gz.fq.gz/}
	echo "Bwa aln, file sort and index for: $output_basename"
        bwa aln -t $NCORES -l $SEED -o $GAP -n $EDIT $REF $input_file > $ALN_OUT_DIR_1/$output_basename.sai
	bwa samse $REF $ALN_OUT_DIR_1/$output_basename.sai $input_file | samtools view -Sb -o $ALN_OUT_DIR_1/$output_basename.bam
        samtools sort -@ $NCORES -m 3G $ALN_OUT_DIR_1/$output_basename.bam -o $ALN_OUT_DIR_1/$output_basename.sorted.bam
        samtools index $ALN_OUT_DIR_1/$output_basename.sorted.bam
done

mkdir -p $ALN_OUT_DIR_2
for input_file in `ls $ALN_OUT_DIR_1/*.sorted.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/.sorted.bam/}
	NAME=$file_basename
	echo "Add read group for: $output_basename"
	java -jar $PICARD_HOME/bin/picard.jar AddOrReplaceReadGroups INPUT=$input_file \
                                            OUTPUT=$ALN_OUT_DIR_2/$output_basename"_RG.bam" \
                                            RGID=$NAME RGLB=$NAME RGPL="Illumina" RGPU="collapsed" RGSM=$NAME VALIDATION_STRINGENCY=LENIENT
        samtools index $ALN_OUT_DIR_2/$output_basename"_RG.bam"
done

mkdir -p $ALN_OUT_DIR_3
mkdir -p $ALN_OUT_DIR_4
for input_file in `ls $ALN_OUT_DIR_2/*_RG.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_RG.bam/}
        NAME=$file_basename
        echo "Removing low quality reads from: $output_basename"
        samtools view -b -h -q 30 $input_file > $ALN_OUT_DIR_3/$file_basename".quality.bam"
        samtools index $ALN_OUT_DIR_3/$file_basename".quality.bam"
        echo "Only primary alignments will be retained for file: $output_basename"
	samtools view -bF 2304 $ALN_OUT_DIR_3/$file_basename".quality.bam" > $ALN_OUT_DIR_4/$file_basename".quality.primary.bam"
        samtools index $ALN_OUT_DIR_4/$file_basename".quality.primary.bam"
done

exit
