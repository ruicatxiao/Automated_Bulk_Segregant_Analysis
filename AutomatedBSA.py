#!/usr/bin/env python3


import argparse
import subprocess
import os
import logging
import sys
from pyfiglet import Figlet
import pandas as pd
from pathlib import Path
from colorama import init, Fore, Style
from tqdm import tqdm


# v2 script, has not tested yet

def setup_logging():
    # Initialize colorama
    init(autoreset=True)
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # Set root logger level to DEBUG

    # Formatter for logs
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')

    # File handler - logs all messages
    file_handler = logging.FileHandler("AutomatedBSA.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler - logs INFO and above
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
        # Stream output to screen and log
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
        # Assuming Trim Galore appends '_val_1.fq.gz' and '_val_2.fq.gz' to the original filenames
        trimmed_r1 = read1.name.replace(".fastq.gz", "_val_1.fq.gz")
        trimmed_r2 = read2.name.replace(".fastq.gz", "_val_2.fq.gz")
        trimmed_read1.append(str(trimmed_out / trimmed_r1))
        trimmed_read2.append(str(trimmed_out / trimmed_r2))
        logging.info(f"Trimmed reads for sample {sample}: {trimmed_r1}, {trimmed_r2}")
    df['trimmed_read1'] = trimmed_read1
    df['trimmed_read2'] = trimmed_read2
    return df


def index_reference(ref_genome, output_dir, pbar):
    logging.info(f"Indexing reference genome: {ref_genome}")
    # bwa index
    run_command(["bwa", "index", str(ref_genome)], pbar=pbar)
    # samtools faidx
    run_command(["samtools", "faidx", str(ref_genome)], pbar=pbar)
    # gatk CreateSequenceDictionary
    run_command([
        "gatk", "CreateSequenceDictionary",
        "-R", str(ref_genome)
    ], pbar=pbar)
    # Increment progress bar for indexing
    # Assuming each command is one task
    # If multiple commands are considered as separate tasks, adjust accordingly
    # Here, already incremented in run_command


def map_reads(df, ref_genome, output_dir, threads, pbar):
    logging.info("Starting read mapping with BWA MEM")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        trimmed_r1 = Path(row['trimmed_read1'])
        trimmed_r2 = Path(row['trimmed_read2'])
        rg = f"@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:{sample}/"
        sam_output = output_dir / 'sam_bam' / f"{sample}.sam"
        command = [
            "bwa", "mem",
            "-t", str(threads),
            str(ref_genome),
            str(trimmed_r1),
            str(trimmed_r2),
            "-M",
            "-R", rg
        ]
        with open(sam_output, 'w') as f:
            logging.info(f"Running command: {' '.join(command)}")
            process = subprocess.Popen(command, stdout=f, stderr=subprocess.PIPE, text=True)
            _, stderr = process.communicate()
            if process.returncode != 0:
                logging.error(f"BWA MEM failed for sample {sample}: {stderr}")
                sys.exit(1)
        logging.info(f"Mapped reads for sample {sample}, SAM saved to {sam_output}")
        if pbar:
            pbar.update(1)


def convert_sort_bam(df, output_dir, threads, pbar):
    logging.info("Converting SAM to BAM and sorting with samtools")
    for idx, row in df.iterrows():
        sample = row['sampleName']
        sam_file = output_dir / 'sam_bam' / f"{sample}.sam"
        bam_file = output_dir / 'sam_bam' / f"{sample}.bam"
        sorted_bam = output_dir / 'sam_bam' / f"{sample}.sorted.bam"
        
        # samtools view
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
        
        # samtools sort
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
    
    # Define filter expressions
    info_filter = "QD > 2.0 & FS < 60.0 & SOR < 3.0"
    genotype_filter = "DP > 10 & GQ > 90"
    
    # Construct the vcffilter command
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
    fai_file = Path(ref_genome).with_suffix('.fasta.fai')
    # Update the path to the custom script in the bin folder
    script_path = Path.cwd() / 'bin' / 'scatter_plot_snp_location.py'
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
    # Update the path to the custom R script in the bin folder
    script_path = Path.cwd() / 'bin' / 'BSA_R_Preprocessing.R'
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
        # Move all .tsv files to the 'tables' folder
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

        # Move all .pdf files to the 'plots' folder
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

    create_directories(output_dir)

    df = read_samplesheet(samplesheet)
    num_samples = len(df)
    logging.info(f"Number of samples to process: {num_samples}")

    # Calculate total tasks: 6 per sample + 16 non-sample steps
    total_tasks = 6 * num_samples + 16

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
        create_vcf_list(output_dir, pbar)
        merge_all_vcfs(ref_genome, create_vcf_list(output_dir, pbar), output_dir, pbar)
        genotype_combined_vcf(ref_genome, output_dir, pbar)
        select_bulk_snps(ref_genome, output_dir, pbar)
        filter_bulk_snps(ref_genome, output_dir, pbar)
        variants_to_table(ref_genome, output_dir, pbar)
        run_post_processing(output_dir, pbar)
        organize_output_files(output_dir, pbar)

    logging.info("Automated BSA Workflow Completed Successfully")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        # No arguments provided
        print_colored_splash()
        print("No arguments provided. Use --help for more information.")
        sys.exit(0)
    main()
