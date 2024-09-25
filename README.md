# Automated Bulk Segregant Analysis (ABSA)

![BSA drawio](https://github.com/user-attachments/assets/ab752971-1ee2-4b3a-97b3-65fa1bbb4173)


## Description

Automated Bulk Segregant Analysis (ABSA) is a pipeline designed to streamline the Bulk Segregant Analysis (BSA) workflow for Cryptosporidium crossings genomic data. It automates a series of bioinformatics tools and processes, enhancing efficiency, reproducibility, and scalability. This pipeline is a streamlined implementation of Dr. Xue Li's BSA analysis. Bulk segregant analysis is a technique for identifying the genetic loci that underlie phenotypic trait differences. The basic approach is to compare two pools of individuals from the opposing tails of the phenotypic distribution, sampled from an interbred population. This pipeline takes input all the way to allele freq table generation and plotting for all valid SNPs.

## Features

- **Read Trimming**: Utilizes Trim Galore for trimming sequencing reads.
- **Read Mapping**: Employs BWA MEM for aligning reads to a reference genome.
- **BAM Processing**: SAM/BAM processing using SAMTools and GATK4.
- **Variant Calling**: Calls variants using GATK4 HaplotypeCaller.
- **SNP Filtering**: Filters SNPs using vcffilter and GATK4 SelectVariants.
- **Post-Processing**: Generates plots and tables for downstream analysis.
- **Organized Outputs**: Automatically organizes output tasble and plot files into designated folders (`tables` and `plots`).
- **Comprehensive Logging**: Logs detailed workflow progress and errors to both the console and a log file (`AutomatedBSA.log`).

## Requirements

### System Dependencies

Ensure the following bioinformatics tools are installed and accessible in your system's `$PATH`:

- [**Trim Galore**](https://github.com/FelixKrueger/TrimGalore)
- [**BWA**](https://github.com/lh3/bwa)
- [**SAMTools**](https://github.com/samtools/samtools)
- [**GATK4**](https://github.com/broadinstitute/gatk)
- [**vcffilter**](https://github.com/vcflib/vcflib)
- [**R**](https://www.r-project.org/)
- **Python 3.10.11**: The pipeline has been tested with Python version 3.10.11.

### Python Dependency

pandas pyfiglet colorama tqdm

### R Dependency

ggplot2 readr

## Installation

- git clone https://github.com/ruicatxiao/Automated_Bulk_Segregant_Analysis.git

- chmod u+x AutomatedBSA.py

- chmod u+x bin/scatter_plot_snp_location.py

- chmod u+x bin/BSA_R_Preprocessing.R

## Prepare data
Follow the provided samplesheet.csv file in the repo and place read files into raw_reads folder, update samplesheet.csv to reflect changes. You do not need to provide absolute path to read files

You should have a reference genome in fasta, a samplesheet.csv and the raw_reads folder containing all reads

## Usage

python3 AutomatedBSA.py \
--ref <GENOME_REFERENCE.fasta> \
--sample samplesheet.csv \
--threads <NUMBER_OF_CPU_THREADS>

- CpBGF genome is provided by default. replace this with any other Cryptosporidum genome as needed


## Output
All final output tables are located in "tables" folder. All generated plots file are located in "plots" folder


## In Development
- Outlier SNPs removal via median absolute deviation (MAD)
- TriCube smoothing of SNP frequency distribution
- QTLSeqR processing, g-prime calculation and plotting
- nf-core adaptation to remove local dependency requirements

## Credits

ABSA is conceptualized by Sebastian Shaw and Rui Xiao. The pipeline is developed and implemented by Rui Xiao.

We thank Xue Li for implementing the original BSA, which ABSA is based upon

## Citations

=====PLACE HOLDER FOR STRIPEN GROUP PUBLICATION CITATION USING ABSA=======

Brenneman KV, Li X, Kumar S, Delgado E, Checkley LA, Shoue DA, Reyes A, Abatiyow BA, Haile MT, Tripura R, Peto T, Lek D, Button-Simons KA, Kappe SHI, Dhorda M, Nosten F, Nkhoma SC, Cheeseman IH, Vaughan AM, Ferdig MT, Anderson TJC. Optimizing bulk segregant analysis of drug resistance using Plasmodium falciparum genetic crosses conducted in humanized mice. iScience. 2022 Mar 16;25(4):104095. doi: 10.1016/j.isci.2022.104095. PMID: 35372813; PMCID: PMC8971943.
