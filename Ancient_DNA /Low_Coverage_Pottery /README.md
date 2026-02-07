# Low Coverage aDNA Analysis

This GitHub repository contains all files and scripts used to analyze samples from the **"Low Coverage"** experiment, focusing on **ancient DNA (aDNA)** data generated from **three archaeological remains**.

---

## Data availability

Raw sequencing data will be deposited in the **NCBI Sequence Read Archive (SRA)** and will be accessible at:  
**[link to be added]**

---

## Repository structure and computing environment

- The Bash scripts used for the bioinformatic analyses are available in the `Bash_Script/` directory.
- All Bash scripts were executed on the **Galileo 100 high-performance computing (HPC) cluster** at **CINECA (Italy)**.
- The R scripts used for plotting and figure generation are available in the `R/` directory.

---

## Software

The following software tools were used (alphabetical order):

| Software      | Version                | Function                                                                 | DOI or LINK |
|---------------|------------------------|---------------------------------------------------------------------------|-----|
| BWA           | 0.7.12                 | Alignment of sequencing reads to the reference genome                     |  10.1093/bioinformatics/btp324   |
| Cutadapt      | 4.1                    | Adapter trimming and removal of low-quality bases                         |  10.14806/ej.17.1.200   |
| FastQC        | 0.11.9                 | Quality control assessment of raw sequencing data                         |  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/   |
| Kraken2       | 2.17.1                 | Taxonomic classification of metagenomic reads                             | https://doi.org/10.1186/s13059-019-1891-0    |
| LeeHom        | 1.2.15                 | Merging of overlapping paired-end reads                                   |  https://doi.org/10.1093/nar/gku699   |
| mapDamage     | 2.2.0-86-g81d0aca      | Estimation and visualization of post-mortem DNA damage patterns           |   10.1093/bioinformatics/btr347  |
| MetaPhlAn     | 4.2.4                  | Profiling of microbial community composition                              | 10.1038/s41587-023-01688-w    |
| Picard        | 2.25.7                 | Manipulation and filtering of SAM/BAM files                               |  http://broadinstitute.github.io/picard/   |
| PMDtools      | 0.50                   | Filtering of ancient DNA reads based on post-mortem damage patterns       | 10.1073/pnas.1318934111    |
| Samtools      | 1.13                   | Processing, indexing, and manipulation of alignment files                 | 10.1093/bioinformatics/btp352   |
| Trimmomatic   | 0.39                   | Quality trimming and filtering of sequencing reads                        |  10.1093/bioinformatics/btu170   |

---

## Reference genomes and databases

### Human reference genome (ENSEMBL)

The human reference genome (**GRCh38**) was downloaded from **ENSEMBL**.

```bash
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
```

### Archaeal reference sequences (NCBI)

Archaeal reference sequences were downloaded from **NCBI**:

- `NZ_LR698975.1`
- `NC_009515.1`

### Kraken2 database

The Kraken2 database was obtained from the official Kraken repository.

```bash
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20241228.tar.gz
```

### MetaPhlAn database

The MetaPhlAn database used in this study was:

- `mpa_vJan25_CHOCOPhlAnSGB_202503.fna`

This database is downloadable from the official MetaPhlAn repository.

# Bioinformatic pipeline

The bioinformatic analyses were performed using a series of custom Bash scripts, described below.


`001_QualityCheck`

This script performs quality control and preprocessing of raw sequencing reads.  
FastQC is used to assess the quality of the raw data. Adapter sequences are removed using Cutadapt.  
Reads are then filtered using Trimmomatic to retain sequences with an average Phred quality score of at least 30 and a minimum length of 30 bp.  
FastQC is run again to evaluate post-filtering quality.  
Finally, overlapping paired-end reads are merged using LeeHom.

`002_Align_Human`

Reads obtained from the `001_QualityCheck` step are aligned to the human reference genome.  
The reference genome is downloaded from Ensembl and indexed using BWA.  
Read alignment is performed using recommended parameters for ancient DNA data.  
Resulting BAM files are sorted and indexed with Samtools, and only high-quality reads are retained.

`003_DNADamagePattern_Human`

Post-mortem DNA damage patterns are evaluated using the BAM files generated in `002_Align_Human`.  
DNA damage is quantified using PMDtools, and damage profiles are visualized using mapDamage.  
Text files containing edit distance information and read counts are generated for each sample.

`004_Align_Bacteria_KrakenBracken`

Reads that did not align to the human reference genome are extracted from the BAM files produced in `002_Align_Human`, and are assumed to be of non-human origin.  
The microbial composition of the samples is explored using Kraken2.

`004_Align_Bacteria_MPA`

Similarly, non-human reads obtained from `002_Align_Human` are used to characterize the microbial profile of the samples using MetaPhlAn4.

`005_Align_Bacteria_P`

Pathogenic taxa identified as relevant and consistently detected by both Kraken2 and MetaPhlAn are selected for further analysis.  
Reference genomes of the selected pathogens are compiled into a FASTA file and indexed using BWA.  
Only reads not aligned to the human genome are mapped against these bacterial reference genomes.  
In this study, the bacterial genomes used were *Methanobrevibacter smithii* and *Methanosphaera stadtmanae*.  
Only high-quality alignments are retained using Samtools.

`006_DNADamagePattern_Bacteria_P`

DNA damage patterns are evaluated on the BAM files generated in `005_Align_Bacteria_P`, following the same strategy described in `003_DNADamagePattern_Human`.  
Text files containing edit distance information and read counts are also generated for each sample.
