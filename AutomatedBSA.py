#!/usr/bin/env python3


import argparse
import subprocess
import os
import logging
import sys
import shutil # Add this import
from pyfiglet import Figlet
import pandas as pd
from pathlib import Path
from colorama import init, Fore, Style
from tqdm import tqdm

# v4 script, debugged, tested working in conda, direct and singularity

def setup_logging():
    init(autoreset=True)
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) 

    # Formatter for logs
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')

    file_handler = logging.FileHandler("AutomatedBSA.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def print_colored_splash():
    f = Figlet(font='slant')
    splash = f.renderText('AUTOMATED BULK SEGREGANT ANALYSIS')
    colored_splash = (
        Fore.CYAN + splash
    )
    print(colored_splash)


def run_command(command, shell=False, pbar=None, update_pbar=True):
    """
    Runs a shell command and logs it. Streams output to the screen and logs.
    Raises an exception if the command fails.
    
    Parameters:
    - command (list): Command and arguments to execute.
    - shell (bool): Whether to execute the command through the shell.
    - pbar (tqdm.tqdm): The progress bar object to update.
    - update_pbar (bool): Whether to update the progress bar after completion.
    """
    logging.info(f"Running command: {' '.join(command) if isinstance(command, list) else command}")
    try:
        process = subprocess.Popen(
            command,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        for line in process.stdout:
            print(line, end='')
            logging.debug(line.strip())
        process.wait()
        if process.returncode != 0:
            logging.error(f"Command failed with return code {process.returncode}")
            sys.exit(1)
        if update_pbar and pbar:
            pbar.update(1)
    except Exception as e:
        logging.error(f"Exception occurred while running command: {e}")
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Automate BSA Workflow",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Path to reference genome fasta"
    )
    parser.add_argument(
        "--sample",
        required=True,
        help="Path to samplesheet.csv"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads to use"
    )
    return parser.parse_args()

def check_required_tools():
    """
    Checks if all required tools are available in the system PATH.
    Exits the script if any tool is missing.
    """
    required_tools = [
        "bwa", "samtools", "gatk", "trim_galore", "vcffilter", "Rscript", "python3"
    ]
    missing_tools = []

    logging.info("Checking for required tools...")
    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)
            logging.error(f"Required tool not found in PATH: {tool}")

    if missing_tools:
        logging.critical(f"Critical Error: The following required tools are missing: {', '.join(missing_tools)}")
        print(f"{Fore.RED}Critical Error: The following required tools are missing: {', '.join(missing_tools)}{Style.RESET_ALL}")
        sys.exit(1)
    else:
        logging.info("All required tools are available in PATH.")


def read_samplesheet(samplesheet_path):
    logging.info(f"Reading samplesheet from {samplesheet_path}")
    df = pd.read_csv(samplesheet_path)
    required_columns = {'sampleName', 'sampleType', 'read1', 'read2'}
    if not required_columns.issubset(df.columns):
        logging.error(f"Samplesheet is missing required columns. Required columns: {required_columns}")
        sys.exit(1)
    return df


def create_directories(output_dir):
    dirs = ['trimmed', 'sam_bam', 'raw_vcf', 'work_vcf', 'tables', 'plots']
    for d in dirs:
        path = output_dir / d
        path.mkdir(parents=True, exist_ok=True)
        logging.info(f"Ensured directory exists: {path}")


def trim_reads(df, output_dir, threads, pbar):
    logging.info("Starting read trimming with Trim Galore")
    trimmed_read1 = []
    trimmed_read2 = []
    for idx, row in df.iterrows():
        sample = row['sampleName']
        read1 = Path(row['read1'])
        read2 = Path(row['read2'])
        trimmed_out = output_dir / 'trimmed'
        
        expected_trimmed_r1_name = read1.name.replace(".fastq.gz", "_val_1.fq.gz")
        expected_trimmed_r2_name = read2.name.replace(".fastq.gz", "_val_2.fq.gz")
        expected_trimmed_r1_path = trimmed_out / expected_trimmed_r1_name
        expected_trimmed_r2_path = trimmed_out / expected_trimmed_r2_name

        if expected_trimmed_r1_path.exists() and expected_trimmed_r1_path.stat().st_size > 0 and \
           expected_trimmed_r2_path.exists() and expected_trimmed_r2_path.stat().st_size > 0:
            logging.info(f"Trimmed files for sample {sample} already exist and are non-empty. Skipping trimming.")
            trimmed_read1.append(str(expected_trimmed_r1_path))
            trimmed_read2.append(str(expected_trimmed_r2_path))
        else:
            command = [
                "trim_galore",
                "--cores", str(threads),
                "--paired",
                "--fastqc",
                "--gzip",
                "-o", str(trimmed_out),
                str(read1),
                str(read2)
            ]
            run_command(command, pbar=pbar)

            trimmed_read1.append(str(expected_trimmed_r1_path))
            trimmed_read2.append(str(expected_trimmed_r2_path))

        logging.info(f"Trimmed reads for sample {sample}: {expected_trimmed_r1_path}, {expected_trimmed_r2_path}")
    df['trimmed_read1'] = trimmed_read1
    df['trimmed_read2'] = trimmed_read2
    return df


def index_reference(ref_genome, output_dir, pbar):
    """
    Indexes the reference genome by running:
    1. bwa index
    2. samtools faidx
    3. gatk CreateSequenceDictionary

    Now includes a safeguard to remove old index files for this reference
    to avoid errors if they already exist.
    """
    logging.info(f"Indexing reference genome: {ref_genome}")

    ref_path = Path(ref_genome).resolve()
    logging.info(f"Resolved reference path: {ref_path}")

    index_extensions_to_remove = [
        '.fai', '.sa', '.amb', '.ann', '.pac', '.bwt', '.dict'
    ]

    for ext in index_extensions_to_remove:

        if ext == '.dict':
            fpath = ref_path.with_suffix(ext)
        else:

            fpath = str(ref_path) + ext
        f = Path(fpath)
        if f.exists():
            logging.info(f"Removing old index file: {f}")
            f.unlink()


    bwa_index_command = ["bwa", "index", str(ref_path)]
    logging.info(f"Running BWA index: {' '.join(bwa_index_command)}")
    run_command(bwa_index_command, pbar=pbar)


    samtools_faidx_command = ["samtools", "faidx", str(ref_path)]
    logging.info(f"Running samtools faidx: {' '.join(samtools_faidx_command)}")
    run_command(samtools_faidx_command, pbar=pbar)

 
    gatk_dict_command = [
        "gatk", "CreateSequenceDictionary",
        "-R", str(ref_path) 
    ]
    logging.info(f"Running GATK CreateSequenceDictionary: {' '.join(gatk_dict_command)}")
    run_command(gatk_dict_command, pbar=pbar)



def map_reads(df, ref_genome, output_dir, threads, pbar):
    """
    Maps trimmed reads to the reference genome using BWA MEM.
    Uses the absolute path of the reference genome file (without extension)
    as the idxbase for BWA.
    """
    logging.info("Starting read mapping with BWA MEM")

    ref_path = Path(ref_genome).resolve() 


    if not ref_path.exists():
        logging.error(f"Reference genome file does not exist: {ref_path}")
        sys.exit(1)


    index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    for ext in index_extensions:
        idx_file = ref_path.with_suffix(ref_path.suffix + ext)
        if not idx_file.exists():
            logging.error(f"Required BWA index file does not exist: {idx_file}")
            logging.error("Please ensure 'bwa index' was run successfully on the reference genome.")
            sys.exit(1)
    logging.info("All required BWA index files found.")

    for idx, row in df.iterrows():
        sample = row['sampleName']
        trimmed_r1 = Path(row['trimmed_read1']).resolve() 
        trimmed_r2 = Path(row['trimmed_read2']).resolve() 
        
        rg = f"@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:{sample}"

        sam_output = output_dir / 'sam_bam' / f"{sample}.sam"


        if sam_output.exists() and sam_output.stat().st_size > 0:
            logging.info(f"SAM file for sample {sample} already exists and is non-empty. Skipping mapping.")
            pbar.update(1) 
            continue 


        command = [
            "bwa", "mem",
            "-t", str(threads),
            "-M", 
            "-R", rg, 
            str(ref_path),  
            str(trimmed_r1), 
            str(trimmed_r2)  
        ]
        
        logging.info(f"Running BWA MEM command for sample {sample}")
        logging.debug(f"Command: {' '.join(command)}") 
        
        try:

            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True #
            )

            with open(sam_output, 'w') as f:
                f.write(result.stdout)
            
            logging.info(f"Mapped reads for sample {sample}, SAM saved to {sam_output}")
            
        except subprocess.CalledProcessError as e:

            error_output = e.stderr if e.stderr else "No stderr output captured."
            logging.error(f"BWA MEM failed for sample {sample} (Return code {e.returncode}):\n{error_output}")
            logging.debug(f"Failed command was: {' '.join(command)}")
            sys.exit(1) 
        except Exception as e:

            logging.error(f"Unexpected exception occurred while running BWA MEM for sample {sample}: {e}")
            logging.debug(f"Failed command was: {' '.join(command)}")
            sys.exit(1)
            
        if pbar:
            pbar.update(1) 


def convert_sort_bam(df, output_dir, threads, pbar):
    logging.info("Converting SAM to BAM and sorting with samtools")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        sam_file = output_dir / 'sam_bam' / f"{sample}.sam"
        bam_file = output_dir / 'sam_bam' / f"{sample}.bam"
        sorted_bam = output_dir / 'sam_bam' / f"{sample}.sorted.bam"
        

        if sorted_bam.exists() and sorted_bam.stat().st_size > 0:
            logging.info(f"Sorted BAM file for sample {sample} already exists and is non-empty. Skipping conversion and sorting.")
            pbar.update(2) 
            continue 
        
        intermediate_bam_exists = bam_file.exists() and bam_file.stat().st_size > 0

        if not intermediate_bam_exists:
            command_view = [
                "samtools", "view",
                "-h",
                "-@", str(threads),
                "-S",
                "-b",
                str(sam_file),
                "-o", str(bam_file)
            ]
            run_command(command_view, pbar=pbar)
            logging.info(f"Converted SAM to BAM for sample {sample}, BAM saved to {bam_file}")
        else:
            logging.info(f"Intermediate BAM file for sample {sample} already exists and is non-empty. Skipping SAM to BAM conversion.")
            pbar.update(1) 

        command_sort = [
            "samtools", "sort",
            "-@", str(threads),
            str(bam_file),
            "-o", str(sorted_bam)
        ]
        run_command(command_sort, pbar=pbar)
        logging.info(f"Sorted BAM for sample {sample} saved to {sorted_bam}")


def mark_duplicates(df, output_dir, pbar):
    logging.info("Marking duplicates with GATK MarkDuplicates")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        sorted_bam = output_dir / 'sam_bam' / f"{sample}.sorted.bam"
        final_bam = output_dir / 'sam_bam' / f"{sample}.final.bam"
        metrics = output_dir / 'sam_bam' / f"{sample}_dup_metrics.txt"
        
        if final_bam.exists() and final_bam.stat().st_size > 0:
            logging.info(f"Final BAM file for sample {sample} already exists and is non-empty. Skipping duplicate marking.")
            pbar.update(1)
            continue 

        command = [
            "gatk", "MarkDuplicates",
            "-I", str(sorted_bam),
            "-O", str(final_bam),
            "-M", str(metrics)
        ]
        run_command(command, pbar=pbar)
        logging.info(f"Marked duplicates for sample {sample}, final BAM: {final_bam}")


def index_bam(df, output_dir, pbar):
    logging.info("Indexing BAM files with samtools index")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        final_bam = output_dir / 'sam_bam' / f"{sample}.final.bam"
        bai_file = final_bam.with_suffix(final_bam.suffix + '.bai') 
        if bai_file.exists() and bai_file.stat().st_size > 0:
            logging.info(f"BAM index file for sample {sample} already exists and is non-empty. Skipping indexing.")
            pbar.update(1)
            continue

        command = [
            "samtools", "index",
            str(final_bam)
        ]
        run_command(command, pbar=pbar)
        logging.info(f"Indexed BAM for sample {sample}")
    

def call_variants(df, ref_genome, output_dir, threads, pbar):
    logging.info("Calling variants with GATK HaplotypeCaller")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        final_bam = output_dir / 'sam_bam' / f"{sample}.final.bam"
        vcf_output = output_dir / 'raw_vcf' / f"{sample}.vcf"
        command = [
            "gatk", "HaplotypeCaller",
            "-R", str(ref_genome),
            "-I", str(final_bam),
            "-ERC", "BP_RESOLUTION",
            "-O", str(vcf_output),
            "-new-qual",
            "-ploidy", "1"
        ]
        run_command(command, pbar=pbar)
        logging.info(f"Called variants for sample {sample}, VCF saved to {vcf_output}")


def identify_and_rename_parents(df, output_dir, pbar):
    logging.info("Identifying parent samples and renaming their VCFs")
    parents = df[df['sampleType'].isin(['parent_1', 'parent_2'])]
    if parents.shape[0] != 2:
        logging.error("Expected exactly two parent samples (parent_1 and parent_2)")
        sys.exit(1)
    for idx, row in parents.iterrows():
        sample = row['sampleName']
        sample_type = row['sampleType']
        original_vcf = output_dir / 'raw_vcf' / f"{sample}.vcf"
        new_vcf = output_dir / 'raw_vcf' / f"{sample_type}.vcf"
        original_vcf.rename(new_vcf)
        logging.info(f"Renamed {original_vcf} to {new_vcf}")
    if pbar:
        pbar.update(1)


def combine_parent_vcfs(ref_genome, output_dir, pbar):
    logging.info("Combining parent VCFs with GATK CombineGVCFs")
    parent_vcf1 = output_dir / 'raw_vcf' / 'parent_1.vcf'
    parent_vcf2 = output_dir / 'raw_vcf' / 'parent_2.vcf'
    combined_vcf = output_dir / 'work_vcf' / 'parents.vcf'
    command = [
        "gatk", "CombineGVCFs",
        "-R", str(ref_genome),
        "-V", str(parent_vcf1),
        "-V", str(parent_vcf2),
        "-O", str(combined_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Combined parent VCFs into {combined_vcf}")


def genotype_parents(ref_genome, output_dir, pbar):
    logging.info("Genotyping parents with GATK GenotypeGVCFs")
    combined_vcf = output_dir / 'work_vcf' / 'parents.vcf'
    genotype_vcf = output_dir / 'work_vcf' / 'p.vcf'
    command = [
        "gatk", "GenotypeGVCFs",
        "-R", str(ref_genome),
        "-V", str(combined_vcf),
        "-O", str(genotype_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Genotyped parents, output VCF: {genotype_vcf}")


def select_parent_snps(ref_genome, output_dir, pbar):
    logging.info("Selecting parental SNPs with GATK SelectVariants")
    genotype_vcf = output_dir / 'work_vcf' / 'p.vcf'
    snp_vcf = output_dir / 'work_vcf' / 'p.SNP.vcf'
    command = [
        "gatk", "SelectVariants",
        "-R", str(ref_genome),
        "-V", str(genotype_vcf),
        "--select-type-to-include", "SNP",
        "-O", str(snp_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Selected SNPs for parents, VCF: {snp_vcf}")
    if pbar:
        pbar.update(1)


def filter_parent_snps(output_dir, pbar):
    logging.info("Filtering parental SNPs with vcffilter")
    input_vcf = output_dir / 'work_vcf' / 'p.SNP.vcf'
    output_vcf = output_dir / 'work_vcf' / 'p.SNP.hardfilter.vcf'
    
    info_filter = "QD > 2.0 & FS < 60.0 & SOR < 3.0"
    genotype_filter = "DP > 10 & GQ > 90"
    
    command = [
        "vcffilter",
        "-f", info_filter,
        "-g", genotype_filter,
        str(input_vcf)
    ]
    
    try:
        with open(output_vcf, 'w') as f:
            logging.info(f"Executing command: {' '.join(command)}")
            subprocess.run(
                command,
                check=True,
                stdout=f,
                stderr=subprocess.PIPE,
                text=True
            )
        logging.info(f"Filtered parental SNPs, output VCF: {output_vcf}")
        if pbar:
            pbar.update(1)
    except subprocess.CalledProcessError as e:
        logging.error(f"vcffilter failed for {input_vcf}: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred during vcffilter: {e}")
        sys.exit(1)


def select_biallelic_variants(ref_genome, output_dir, pbar):
    logging.info("Selecting biallelic variants with GATK SelectVariants")
    input_vcf = output_dir / 'work_vcf' / 'p.SNP.hardfilter.vcf'
    output_vcf = output_dir / 'work_vcf' / 'p.SNP.valid.vcf'
    command = [
        "gatk", "SelectVariants",
        "-R", str(ref_genome),
        "-V", str(input_vcf),
        "--exclude-non-variants", "true",
        "--max-nocall-number", "0",
        "--remove-unused-alternates", "true",
        "--restrict-alleles-to", "BIALLELIC",
        "--select", "AF > 0.25 && AF < 1.0",
        "-O", str(output_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Selected biallelic variants, output VCF: {output_vcf}")
    if pbar:
        pbar.update(1)



def plot_snp_distribution(output_dir, ref_genome, pbar):
    logging.info("Plotting SNP distribution with scatter_plot_snp_location.py")
    vcf_file = output_dir / 'work_vcf' / 'p.SNP.valid.vcf'
    fai_file = Path(str(ref_genome) + ".fai")

    script_path = Path(__file__).parent / 'scatter_plot_snp_location.py'
    
    command = [
        "python3", str(script_path),
        str(vcf_file),
        str(fai_file)
    ]
    run_command(command, pbar=pbar)
    logging.info("Completed SNP distribution plot")
    if pbar:
        pbar.update(1)


def create_vcf_list(output_dir, pbar):
    logging.info("Creating list of all VCF files for merging")
    raw_vcf_dir = output_dir / 'raw_vcf'
    list_file = raw_vcf_dir / 'gvcf.list'
    vcf_files = list(raw_vcf_dir.glob("*.vcf"))
    with open(list_file, 'w') as f:
        for vcf in vcf_files:
            f.write(f"{vcf}\n")
    logging.info(f"Created VCF list file at {list_file}")
    if pbar:
        pbar.update(1)
    return list_file


def merge_all_vcfs(ref_genome, vcf_list, output_dir, pbar):
    logging.info("Merging all VCFs with GATK CombineGVCFs")
    combined_vcf = output_dir / 'work_vcf' / 'combined.vcf'
    command = [
        "gatk", "CombineGVCFs",
        "-R", str(ref_genome),
        "-V", str(vcf_list),
        "-O", str(combined_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Merged VCFs into {combined_vcf}")
    if pbar:
        pbar.update(1)


def genotype_combined_vcf(ref_genome, output_dir, pbar):
    logging.info("Genotyping combined VCF with GATK GenotypeGVCFs")
    combined_vcf = output_dir / 'work_vcf' / 'combined.vcf'
    genotype_vcf = output_dir / 'work_vcf' / 'combined_genotyped.vcf'
    command = [
        "gatk", "GenotypeGVCFs",
        "-R", str(ref_genome),
        "-V", str(combined_vcf),
        "-O", str(genotype_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Genotyped combined VCF, output: {genotype_vcf}")
    if pbar:
        pbar.update(1)


def select_bulk_snps(ref_genome, output_dir, pbar):
    logging.info("Selecting bulk SNPs with GATK SelectVariants")
    combined_genotyped_vcf = output_dir / 'work_vcf' / 'combined_genotyped.vcf'
    bulk_snp_vcf = output_dir / 'work_vcf' / 'combined.SNP.vcf'
    command = [
        "gatk", "SelectVariants",
        "-R", str(ref_genome),
        "-V", str(combined_genotyped_vcf),
        "--select-type-to-include", "SNP",
        "-O", str(bulk_snp_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Selected bulk SNPs, output VCF: {bulk_snp_vcf}")
    if pbar:
        pbar.update(1)


def filter_bulk_snps(ref_genome, output_dir, pbar):
    logging.info("Filtering bulk SNPs with parental SNPs")
    bulk_snp_vcf = output_dir / 'work_vcf' / 'combined.SNP.vcf'
    parent_valid_vcf = output_dir / 'work_vcf' / 'p.SNP.valid.vcf'
    filtered_vcf = output_dir / 'work_vcf' / 'combined.SNP.filtered.vcf'
    command = [
        "gatk", "SelectVariants",
        "-R", str(ref_genome),
        "-V", str(bulk_snp_vcf),
        "--concordance", str(parent_valid_vcf),
        "--restrict-alleles-to", "BIALLELIC",
        "-O", str(filtered_vcf)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Filtered bulk SNPs, output VCF: {filtered_vcf}")
    if pbar:
        pbar.update(1)


def variants_to_table(ref_genome, output_dir, pbar):
    logging.info("Converting variants to table with GATK VariantsToTable")
    filtered_vcf = output_dir / 'work_vcf' / 'combined.SNP.filtered.vcf'
    table_output = output_dir / 'tables' / 'FINAL_SNP_filtered.tsv'
    command = [
        "gatk", "VariantsToTable",
        "-R", str(ref_genome),
        "-V", str(filtered_vcf),
        "-F", "CHROM",
        "-F", "POS",
        "-F", "REF",
        "-F", "ALT",
        "-GF", "AD",
        "-GF", "DP",
        "-GF", "GQ",
        "-GF", "PL",
        "-O", str(table_output)
    ]
    run_command(command, pbar=pbar)
    logging.info(f"Converted variants to table, output: {table_output}")
    if pbar:
        pbar.update(1)



def run_post_processing(output_dir, pbar):
    logging.info("Running post-processing with BSA_R_Preprocessing.R")
    table_file = output_dir / 'tables' / 'FINAL_SNP_filtered.tsv'

    script_path = Path(__file__).parent / 'BSA_R_Preprocessing.R'
    
    command = [
        "Rscript",
        str(script_path),
        str(table_file)
    ]
    run_command(command, pbar=pbar)
    logging.info("Completed post-processing with BSA_R_Preprocessing.R")
    if pbar:
        pbar.update(1)


def organize_output_files(output_dir, pbar):
    logging.info("Organizing output files into 'tables' and 'plots' folders")

    tables_dir = output_dir / 'tables'
    plots_dir = output_dir / 'plots'

    try:

        tsv_files = list(output_dir.rglob("*.tsv"))
        if tsv_files:
            for tsv in tsv_files:
                # Skip if the file is already in the 'tables' directory
                if tsv.parent != tables_dir:
                    destination = tables_dir / tsv.name
                    tsv.rename(destination)
                    logging.info(f"Moved TSV file: {tsv} -> {destination}")
        else:
            logging.info("No .tsv files found to move.")


        pdf_files = list(output_dir.rglob("*.pdf"))
        if pdf_files:
            for pdf in pdf_files:
                # Skip if the file is already in the 'plots' directory
                if pdf.parent != plots_dir:
                    destination = plots_dir / pdf.name
                    pdf.rename(destination)
                    logging.info(f"Moved PDF file: {pdf} -> {destination}")
        else:
            logging.info("No .pdf files found to move.")

        logging.info("Successfully organized output files into 'tables' and 'plots' folders.")
        if pbar:
            pbar.update(1)
    except Exception as e:
        logging.error(f"An error occurred while organizing output files: {e}")
        sys.exit(1)

def print_colored_goodbye():
    f = Figlet(font='slant')
    goodbye_text = 'ABSA COMPLETE - GOODBYE!'
    splash = f.renderText(goodbye_text)
    colored_splash = Fore.GREEN + splash
    print(colored_splash)

def main():
    setup_logging()
    print_colored_splash()
    args = parse_arguments()
    ref_genome = Path(args.ref).resolve()
    samplesheet = Path(args.sample).resolve()
    THREADS = args.threads

    logging.info(f"Reference Genome: {ref_genome}")
    logging.info(f"Samplesheet: {samplesheet}")
    logging.info(f"CPU Threads: {THREADS}")

    output_dir = Path.cwd()

    if not ref_genome.exists():
        logging.error(f"Reference genome not found at {ref_genome}")
        sys.exit(1)
    if not samplesheet.exists():
        logging.error(f"Samplesheet not found at {samplesheet}")
        sys.exit(1)

    # Check for required tools before proceeding
    check_required_tools()
    create_directories(output_dir)
    df = read_samplesheet(samplesheet)
    num_samples = len(df)
    logging.info(f"Number of samples to process: {num_samples}")

    # Calculate total tasks: 6 per sample + 18 non-sample steps (index_reference is now 3 steps counted separately)
    # Original calculation was 6 per sample + 16 fixed steps. Indexing now counts 3 steps in run_command.
    # So it's 6 per sample + (16 - 1 + 3) = 6 * num_samples + 18
    total_tasks = 6 * num_samples + 18

    with tqdm(total=total_tasks, desc='Automated BSA Progress', unit='task') as pbar:
        df = trim_reads(df, output_dir, THREADS, pbar)
        index_reference(ref_genome, output_dir, pbar) 
        map_reads(df, ref_genome, output_dir, THREADS, pbar) 
        convert_sort_bam(df, output_dir, THREADS, pbar)
        mark_duplicates(df, output_dir, pbar)
        index_bam(df, output_dir, pbar)
        call_variants(df, ref_genome, output_dir, THREADS, pbar)
        identify_and_rename_parents(df, output_dir, pbar)
        combine_parent_vcfs(ref_genome, output_dir, pbar)
        genotype_parents(ref_genome, output_dir, pbar)
        select_parent_snps(ref_genome, output_dir, pbar)
        filter_parent_snps(output_dir, pbar)
        select_biallelic_variants(ref_genome, output_dir, pbar)
        plot_snp_distribution(output_dir, ref_genome, pbar)
        vcf_list_path = create_vcf_list(output_dir, pbar) 
        merge_all_vcfs(ref_genome, vcf_list_path, output_dir, pbar) 
        genotype_combined_vcf(ref_genome, output_dir, pbar)
        select_bulk_snps(ref_genome, output_dir, pbar)
        filter_bulk_snps(ref_genome, output_dir, pbar)
        variants_to_table(ref_genome, output_dir, pbar)
        run_post_processing(output_dir, pbar)
        organize_output_files(output_dir, pbar)

    print_colored_goodbye() 
    print(f"{Fore.GREEN}ABSA has completed successfully, Goodbye!{Style.RESET_ALL}")
    logging.info("Successfully completed full pipeline ")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        # No arguments provided
        print_colored_splash()
        print("No arguments provided. Use --help for more information.")
        sys.exit(0)
    main()
