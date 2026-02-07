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

module load profile/advanced
module load profile/bioinf
module load autoload samtools/1.13
module load autoload illumina_genome_Homo_sapiens/hg38

################# PATH & VARIABILES ################################################

REF="004_DB_HUMAN/genome.fa"
ALN_OUT_DIR_4="005_ALN_OUT_DIR_HUMAN/003_PRIMARY/"
PMDTOOLS_1="006_PMDTools_HUMAN/001_PMD_1"
PMDTOOLS_2="006_PMDTools_HUMAN/002_PMD_2"
PMDTOOLS_3="006_PMDTools_HUMAN/003_PMD_3"
MAPDAMAGE_1="007_MapDamage_HUMAN/000_PMD_1"
MAPDAMAGE_2="007_MapDamage_HUMAN/000_PMD_2"
MAPDAMAGE_3="007_MapDamage_HUMAN/000_PMD_3"

####################################################################################

source activate /../../../../miniconda3/envs/py27

mkdir -p $PMDTOOLS_1
mkdir -p $PMDTOOLS_2
mkdir -p $PMDTOOLS_3
for input_file in `ls $ALN_OUT_DIR_4/*_RG.bam.quality.primary.bam`; do
  file_basename=$(basename $input_file)
  output_basename=${file_basename/_RG.bam.quality.primary.bam/}
  echo "$file_basename"
  echo "file: $output_basename"
	echo "PMD=1"
  samtools view -h $input_file | pmdtools --threshold 1 --header > $PMDTOOLS_1/temp.sam
  samtools view -Sb $PMDTOOLS_1/temp.sam > $PMDTOOLS_1/$output_basename"_PMD_1.bam"
	rm -r $PMDTOOLS_1/temp.sam
  echo "PMD=2"
	samtools view -h $input_file | pmdtools --threshold 2 --header > $PMDTOOLS_2/temp.sam
  samtools view -Sb $PMDTOOLS_2/temp.sam > $PMDTOOLS_2/$output_basename"_PMD_2.bam"
	rm -r $PMDTOOLS_2/temp.sam
  echo "PMD=3"
	samtools view -h $input_file | pmdtools --threshold 3 --header > $PMDTOOLS_3/temp.sam
  samtools view -Sb $PMDTOOLS_3/temp.sam > $PMDTOOLS_3/$output_basename"_PMD_3.bam"
	rm -r $PMDTOOLS_3/temp.sam
done

conda deactivate


source activate /../../../../miniconda3/envs/mapdamage

mkdir -p $MAPDAMAGE_1
mkdir -p $MAPDAMAGE_2
mkdir -p $MAPDAMAGE_3

# Mapdamage on PMD=1
for input_file in `ls $PMDTOOLS_1/*_PMD_1.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_PMD_1.bam/}
        echo "$file_basename"
        echo "file: $output_basename"
        echo "PMD=1"
	mapDamage -i $input_file -r $REF -d $output_basename --merge-reference-sequences
	mv $output_basename $MAPDAMAGE_1
done

# Mapdamage on PMD=2
for input_file in `ls $PMDTOOLS_2/*_PMD_2.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_PMD_2.bam/}
        echo "$file_basename"                                                                                                                                      #echo "file: $output_basename"
        echo "PMD=2"
        mapDamage -i $input_file -r $REF -d $output_basename --merge-reference-sequences
        mv $output_basename $MAPDAMAGE_2
done

# Mapdamage on PMD=3
for input_file in `ls $PMDTOOLS_3/*_PMD_3.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_PMD_3.bam/}
        echo "$file_basename"                                                                                                                                      #echo "file: $output_basename"
        echo "PMD=3"
        mapDamage -i $input_file -r $REF -d $output_basename --merge-reference-sequences
        mv $output_basename $MAPDAMAGE_3
done

conda deactivate

##### EVALUATE EDIT DISTANCE

echo -e "filename\tPMD\tNumberRead\tEditDistance" > edit_distance_HUMAN.txt

for input_file in `ls $ALN_OUT_DIR_4/*_RG.bam.quality.primary.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_RG.bam.quality.primary.bam/}

	### input files
        file_pmd_0=$input_file
        file_pmd_1=$PMDTOOLS_1/$output_basename"_PMD_1.bam"
        file_pmd_2=$PMDTOOLS_2/$output_basename"_PMD_2.bam"
        file_pmd_3=$PMDTOOLS_3/$output_basename"_PMD_3.bam"

	### compute edit distance, PMD = 0
        samtools view $file_pmd_0 | awk '{for(i=1;i<=NF;i++){print $i}}' | grep "NM" | sed 's/NM:i://g' | sort | uniq -c | 
		awk -v fname="$output_basename" '{print fname "\t0\t" $1 "\t" $2}'  >> edit_distance_HUMAN.txt

	### compute edit distance, PMD = 1
        samtools view $file_pmd_1 | awk '{for(i=1;i<=NF;i++){print $i}}' | grep "NM" | sed 's/NM:i://g' | sort | uniq -c |
		awk -v fname="$output_basename" '{print fname "\t1\t" $1 "\t" $2}'  >> edit_distance_HUMAN.txt

	### compute edit distance, PMD = 2
        samtools view $file_pmd_2 | awk '{for(i=1;i<=NF;i++){print $i}}' | grep "NM" | sed 's/NM:i://g' | sort | uniq -c |
		awk -v fname="$output_basename" '{print fname "\t2\t" $1 "\t" $2}'  >> edit_distance_HUMAN.txt

        ### compute edit distance, PMD = 3
        samtools view $file_pmd_3 | awk '{for(i=1;i<=NF;i++){print $i}}' | grep "NM" | sed 's/NM:i://g' | sort | uniq -c |
		awk -v fname="$output_basename" '{print fname "\t3\t" $1 "\t" $2}'  >> edit_distance_HUMAN.txt
done

### EVALUATE THE NUMBER OF SEQUENCES ALIGNED FOR EACH SAMPLE AND AT EACH STEP

echo -e "filename\tPMD\tNumberRead" > number_of_reads_HUMAN.txt

for input_file in `ls $ALN_OUT_DIR_4/*_RG.bam.quality.primary.bam`; do
        file_basename=$(basename $input_file)
        output_basename=${file_basename/_RG.bam.quality.primary.bam/}

        ### input files
        file_pmd_0=$input_file
        file_pmd_1=$PMDTOOLS_1/$output_basename"_PMD_1.bam"
        file_pmd_2=$PMDTOOLS_2/$output_basename"_PMD_2.bam"
        file_pmd_3=$PMDTOOLS_3/$output_basename"_PMD_3.bam"

        ### compute the number of aligned reads, PMD=0
	read_pmd_0=$(samtools idxstats $file_pmd_0 | awk '{sum+=$3} END {print sum}')
	echo -e "$output_basename\t0\t$read_pmd_0" >> number_of_reads_HUMAN.txt

	### compute the number of aligned reads, PMD=1
        read_pmd_1=$(samtools idxstats $file_pmd_1 | awk '{sum+=$3} END {print sum}')
        echo -e "$output_basename\t1\t$read_pmd_1" >> number_of_reads_HUMAN.txt

	### compute the number of aligned reads, PMD=2
        read_pmd_2=$(samtools idxstats $file_pmd_2 | awk '{sum+=$3} END {print sum}')
        echo -e "$output_basename\t2\t$read_pmd_2" >> number_of_reads_HUMAN.txt

        ### compute the number of aligned reads, PMD=3
	read_pmd_3=$(samtools idxstats $file_pmd_3 | awk '{sum+=$3} END {print sum}')
        echo -e "$output_basename\t3\t$read_pmd_3" >> number_of_reads_HUMAN.txt
done
