# Automated Bulk Segregant Analysis (ABSA)

![Image](https://github.com/user-attachments/assets/ee72028a-5919-40c3-a654-0b8e1aa14600)

## Description

Automated Bulk Segregant Analysis (ABSA) is a pipeline designed to streamline the Bulk Segregant Analysis (BSA) workflow for Cryptosporidium crossings genomic data. It automates a series of bioinformatics tools and processes, enhancing efficiency, reproducibility, and scalability. This pipeline is a streamlined implementation of Dr. Xue Li's BSA analysis. Bulk segregant analysis is a technique for identifying the genetic loci that underlie phenotypic trait differences. The basic approach is to compare two pools of individuals from the opposing tails of the phenotypic distribution, sampled from an interbred population. This pipeline takes input all the way to allele freq table generation and plotting for all valid SNPs.

Version update Jan-30-2025: now added singularity installation instructions, which takes care of all depdendency issues by constructing singulatity image. 

```text
    ___   __  ____________  __  ______  ________________ 
   /   | / / / /_  __/ __ \/  |/  /   |/_  __/ ____/ __ \
  / /| |/ / / / / / / / / / /|_/ / /| | / / / __/ / / / /
 / ___ / /_/ / / / / /_/ / /  / / ___ |/ / / /___/ /_/ / 
/_/  |_\____/ /_/  \____/_/  /_/_/  |_/_/ /_____/_____/  
                                                         
    ____  __  ____    __ __
   / __ )/ / / / /   / //_/
  / __  / / / / /   / ,<   
 / /_/ / /_/ / /___/ /| |  
/_____/\____/_____/_/ |_|  
                           
   _____ ________________  _______________    _   ________
  / ___// ____/ ____/ __ \/ ____/ ____/   |  / | / /_  __/
  \__ \/ __/ / / __/ /_/ / __/ / / __/ /| | /  |/ / / /   
 ___/ / /___/ /_/ / _, _/ /___/ /_/ / ___ |/ /|  / / /    
/____/_____/\____/_/ |_/_____/\____/_/  |_/_/ |_/ /_/     
                                                          
    ___    _   _____    ____  _______ _________
   /   |  / | / /   |  / /\ \/ / ___//  _/ ___/
  / /| | /  |/ / /| | / /  \  /\__ \ / / \__ \ 
 / ___ |/ /|  / ___ |/ /___/ /___/ // / ___/ / 
/_/  |_/_/ |_/_/  |_/_____/_//____/___//____/  
```

## Features

- **Read Trimming**: Utilizes Trim Galore for trimming sequencing reads.
- **Read Mapping**: Employs BWA MEM for aligning reads to a reference genome.
- **BAM Processing**: SAM/BAM processing using SAMTools and GATK4.
- **Variant Calling**: Calls variants using GATK4 HaplotypeCaller.
- **SNP Filtering**: Filters SNPs using vcffilter and GATK4 SelectVariants.
- **Post-Processing**: Generates plots and tables for downstream analysis.
- **Organized Outputs**: Automatically organizes output tasble and plot files into designated folders (`tables` and `plots`).
- **Comprehensive Logging**: Logs detailed workflow progress and errors to both the console and a log file (`AutomatedBSA.log`).

## Installations and Requirements

### Manual Installation
Good for fine tuning specific parameters inside the AutomatedBSA.py

#### System Dependencies

Ensure the following bioinformatics tools are installed and accessible in your system's `$PATH`:

- [**Trim Galore**](https://github.com/FelixKrueger/TrimGalore)
- [**BWA**](https://github.com/lh3/bwa)
- [**SAMTools**](https://github.com/samtools/samtools)
- [**GATK4**](https://github.com/broadinstitute/gatk)
- [**vcffilter**](https://github.com/vcflib/vcflib)
- [**R**](https://www.r-project.org/)
- **Python 3.10.11**: The pipeline has been tested with Python version 3.10.11.

#### Python Dependency

pandas pyfiglet colorama tqdm matplotlib

#### R Dependency

ggplot2 readr

#### Manual Installation Steps
```bash
git clone https://github.com/ruicatxiao/Automated_Bulk_Segregant_Analysis.git
```
```bash
cd Automated_Bulk_Segregant_Analysis/
```
```bash
chmod u+x AutomatedBSA.py
```
```bash
chmod u+x scatter_plot_snp_location.py
```
```bash
chmod u+x BSA_R_Preprocessing.R
```

### Singularity Installation
Good for rapid deployment in new systems

#### System Dependencies
Make sure singularity (>=v3.0.0) is installed on your system and in your $PATH

- [**Singularity**](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

```bash
sudo apt-get install \
   autoconf \
   automake \
   cryptsetup \
   fuse \
   fuse2fs \
   git \
   libfuse-dev \
   libglib2.0-dev \
   libseccomp-dev \
   libtool \
   pkg-config \
   runc \
   squashfs-tools \
   squashfs-tools-ng \
   uidmap \
   wget \
   zlib1g-dev
```

```bash
wget https://github.com/sylabs/singularity/releases/download/v4.2.1/singularity-ce_4.2.1-jammy_amd64.deb
```

```bash
sudo dpkg --install singularity-ce_4.2.1-jammy_amd64.deb

sudo apt-get install -yf

```



#### Automatic Installation Steps
```bash
git clone https://github.com/ruicatxiao/Automated_Bulk_Segregant_Analysis.git
```
```bash
cd Automated_Bulk_Segregant_Analysis/
```
```bash
chmod u+x AutomatedBSA.py
```
```bash
chmod u+x scatter_plot_snp_location.py
```
```bash
chmod u+x BSA_R_Preprocessing.R
```
```bash
sudo singularity build absa.sif absa.def
```

 - (OPTIONAL) Consider adding your singularity image absa.sif to $PATH so it can be executed from everywhere

```bash
mkdir -p ~/singularity/images/
```
```bash
mv absa.sif ~/singularity/images/
```
```bash
echo 'export PATH="$PATH:$HOME/singularity/images/"' >> ~/.bashrc
```
```bash
source ~/.bashrc
```

### Conda Installation
Good for desktop  / laptop installation that requires easier package version controls

#### System Dependencies
Make sure miniconda is installed

#### Create conda env for ABSA
```bash
conda create -n absa python=3.10.11

conda activate absa


conda install \
    trim-galore \
    bwa \
    samtools \
    gatk4 \
    vcflib \
    r-base

conda install \
    pandas \
    matplotlib

pip install pyfiglet colorama tqdm


conda install \
    r-ggplot2 \
    r-readr

```

#### Clone and set up ABSA
```bash
git clone https://github.com/ruicatxiao/Automated_Bulk_Segregant_Analysis.git

cd Automated_Bulk_Segregant_Analysis/

chmod u+x AutomatedBSA.py
chmod u+x scatter_plot_snp_location.py
chmod u+x BSA_R_Preprocessing.R

```

#### Verify installations
```bash
which python
python --version

trim_galore --version
bwa
samtools --version
gatk --version
vcffilter -h
R --version

```

## Usage
- T2TCpBGF genome is provided by default. replace this with any other Cryptosporidum genome as needed

### Manual and Conda Installation Execution
```bash
python3 AutomatedBSA.py \
    --ref <GENOME_REFERENCE.fasta> \
    --sample samplesheet.csv \
    --threads <NUMBER_OF_CPU_THREADS> \
 > /dev/null 2>&1 & # Redirect output and run in background
```

### Manual and Conda Installation Execution in the background
```bash
nohup python3 AutomatedBSA.py \
    --ref <GENOME_REFERENCE.fasta> \
    --sample samplesheet.csv \
    --threads <NUMBER_OF_CPU_THREADS>
```

### Singularity Image Execution
```bash
singularity exec \
    --bind <INPUT_OUTPUT_FOLDER>:/input \       # Mount host directory to container
    absa.sif \                                  # Singularity image
    /bin/bash -c "cd /input && python3 /opt/ABSA/AutomatedBSA.py \
      --ref <GENOME_REFERENCE.fasta> \
      --sample samplesheet.csv \
      --threads <NUMBER_OF_CPU_THREADS>"
```

### Singularity Image Execution in the background
```bash
nohup singularity exec \
    --bind <INPUT_OUTPUT_FOLDER>:/input \        # Mount host directory to container
    absa.sif \                                   # Singularity image
    /bin/bash -c "cd /input && python3 /opt/ABSA/AutomatedBSA.py \
      --ref <GENOME_REFERENCE.fasta> \
      --sample samplesheet.csv \
      --threads <NUMBER_OF_CPU_THREADS>" > /dev/null 2>&1 &       # Redirect output and run in background
```


## Input Structure
Follow the provided samplesheet.csv file in the repo and place read files into raw_reads folder, update samplesheet.csv to reflect changes. You do not need to provide absolute path to read files

You should have a reference genome in fasta, a samplesheet.csv and the raw_reads folder containing all reads

- Manual and Conda Installation Required Input Folder Structure

```text
Automated_Bulk_Segregant_Analysis/
├── AutomatedBSA.py                         # Main workflow script
├── scatter_plot_snp_location.py            # Plotting script
├── BSA_R_Preprocessing.R                   # R Preprocessing script
├── samplesheet.csv                         # Sample sheet
├── refgenome.fasta                         # Reference genome
└── raw_reads/                              # Folder with raw reads
    ├── p1_read1.fastq.gz
    ├── p1_read2.fastq.gz
    ├── p2_read1.fastq.gz
    ├── p2_read2.fastq.gz
    ├── c1_read1.fastq.gz
    ├── c1_read2.fastq.gz
    ├── c2_read1.fastq.gz
    ├── c2_read2.fastq.gz
    ├── c3_read1.fastq.gz
    ├── c3_read2.fastq.gz
    └── ...
```

- Singularity Installation Required Input Folder Structure

```text
INPUT_OUTPUT_FOLDER/
├── samplesheet.csv                         # Sample sheet
├── refgenome.fasta                         # Reference genome
└── raw_reads/                              # Folder with raw reads
    ├── p1_read1.fastq.gz
    ├── p1_read2.fastq.gz
    ├── p2_read1.fastq.gz
    ├── p2_read2.fastq.gz
    ├── c1_read1.fastq.gz
    ├── c1_read2.fastq.gz
    ├── c2_read1.fastq.gz
    ├── c2_read2.fastq.gz
    ├── c3_read1.fastq.gz
    ├── c3_read2.fastq.gz
    └── ...
```

- samplesheet.csv structure
Follow the format below to prepare your samplesheet.csv, as well as reads path
```text

sampleName,sampleType,read1,read2
Parent_1_Name,parent_1,raw_reads/p1_read1.fastq.gz,raw_reads/p1_read2.fastq.gz
Parent_2_Name,parent_2,raw_reads/p2_read1.fastq.gz,raw_reads/p2_read2.fastq.gz
Crossing_Gen1,bulk,raw_reads/c1_read1.fastq.gz,raw_reads/c1_read2.fastq.gz
Crossing_Gen2,bulk,raw_reads/c2_read1.fastq.gz,raw_reads/c2_read2.fastq.gz
Crossing_Gen3,bulk,raw_reads/c3_read1.fastq.gz,raw_reads/c3_read2.fastq.gz
Crossing_Gen4,bulk,raw_reads/c4_read1.fastq.gz,raw_reads/c4_read2.fastq.gz

```

- manuscript_scripts
These are scripts used for generating the figures of the primary manuscript


## Output Structure
All final output tables are located in "tables" folder. All generated plots file are located in "plots" folder

```text
├── plots/              # Final visualizations
├── tables/             # Processed data tables
├── sam_bam/            # Alignment files
├── raw_vcf/            # Initial variant calls
├── work_vcf/           # Filtered variants
└── AutomatedBSA.log    # Process log
```


## In Development
- Adding more functions to R post-processing scripts


## Abstract
The relationship between parasite genotype and pathogenesis is largely unknown for Cryptosporidium, a leading cause of diarrheal disease in children. An array of parasites with similar genomes produces varied disease outcomes in different hosts. Here, we isolate and characterize Cryptosporidium parvum strains that show marked differences in virulence and persistence in mice. Taking advantage of the sexual lifecycle of this eukaryotic pathogen, we use genetic crosses to discover the underlying chromosomal loci. Whole-genome sequencing and bulk segregant analysis of infection selected progeny mapped three loci on chromosomes 2, 6, and 7 associated with the ability to colonize and persist in mice and the positions of drug resistance genes. The chromosome 6 locus encodes the hyper-polymorphic surface glycoprotein GP60. Reverse genetic studies in both parental strains demonstrate that GP60 controls parasite burden and virulence, but not persistence, and reveal the dominance of the less virulent allele, suggesting it restricts virulence.

## Credits

ABSA is conceptualized by Sebastian Shaw and Rui Xiao. The pipeline is developed and implemented by Rui Xiao.

Original BSA algorithm was developped by Xue Li.

## Citations

Shaw, S., Li, X., Buenconsejo, G. Y., Tiffany, Z. H., Cohen, A., Yasur, D., Xiao, R., Beiting, D. P., Anderson, T. J., & Striepen, B. (2025). Genetic crosses reveal genomic loci responsible for virulence in Cryptosporidium parvum infection. bioRxiv (Cold Spring Harbor Laboratory). https://doi.org/10.1101/2025.05.20.655157

Brenneman KV, Li X, Kumar S, Delgado E, Checkley LA, Shoue DA, Reyes A, Abatiyow BA, Haile MT, Tripura R, Peto T, Lek D, Button-Simons KA, Kappe SHI, Dhorda M, Nosten F, Nkhoma SC, Cheeseman IH, Vaughan AM, Ferdig MT, Anderson TJC. Optimizing bulk segregant analysis of drug resistance using Plasmodium falciparum genetic crosses conducted in humanized mice. iScience. 2022 Mar 16;25(4):104095. doi: 10.1016/j.isci.2022.104095. PMID: 35372813; PMCID: PMC8971943.
