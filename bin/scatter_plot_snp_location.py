#!/usr/bin/env python3
"""
Plot SNP Coordinates from a VCF File Across Chromosomes.

This script reads SNP positions from a VCF file and chromosome lengths from a .fai file,
then plots the SNPs across chromosomes scaled to their actual lengths. The output is a PDF
file named 'parents_SNP_filtered_location_plot.pdf'.

Usage:
    python plot_snps.py parental_snp_valid.vcf genome.fasta.fai

Dependencies:
    - Python 3
    - matplotlib

Install Dependencies:
    pip install matplotlib
"""

import matplotlib.pyplot as plt
import sys
import os

def main(vcf_file, fai_file):
    # Define the fixed output file name
    output_file = "parents_SNP_filtered_location_plot.pdf"

    # Verify input files exist
    if not os.path.isfile(vcf_file):
        print(f"Error: VCF file '{vcf_file}' does not exist.")
        sys.exit(1)
    if not os.path.isfile(fai_file):
        print(f"Error: FAI file '{fai_file}' does not exist.")
        sys.exit(1)

    # Read chromosome lengths from the .fai file into a dictionary
    chrom_lengths = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                chrom_id = fields[0]
                try:
                    chrom_length = int(fields[1])
                    chrom_lengths[chrom_id] = chrom_length
                except ValueError:
                    print(f"Warning: Invalid chromosome length for '{chrom_id}' in .fai file.")
    
    # Read SNP coordinates from the VCF file and group them by chromosome
    chrom_coords = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue  # Skip header lines
            fields = line.split('\t')
            if len(fields) >= 2:
                chrom = fields[0]
                try:
                    position = int(fields[1])
                except ValueError:
                    print(f"Warning: Invalid SNP position '{fields[1]}' on chromosome '{chrom}'. Skipping.")
                    continue
                # Include all chromosomes present in the VCF
                if chrom not in chrom_coords:
                    chrom_coords[chrom] = []
                chrom_coords[chrom].append(position)
    
    # Check if any SNPs were found
    if not chrom_coords:
        print("Error: No SNPs found in the VCF file.")
        sys.exit(1)
    
    # Ensure all chromosomes in VCF have corresponding lengths in .fai
    missing_chroms = [chrom for chrom in chrom_coords if chrom not in chrom_lengths]
    if missing_chroms:
        print(f"Warning: Chromosome lengths not found for chromosomes: {', '.join(missing_chroms)}. These chromosomes will be skipped.")
        for chrom in missing_chroms:
            del chrom_coords[chrom]
    
    # If no chromosomes remain after filtering, exit
    if not chrom_coords:
        print("Error: No SNPs to plot after filtering chromosomes without lengths.")
        sys.exit(1)
    
    # Get unique chromosome IDs in sorted order (optional: customize sorting as needed)
    unique_chroms = sorted(chrom_coords.keys())
    
    # Map chromosome IDs to numerical y-axis positions
    chrom_to_y = {chrom: idx for idx, chrom in enumerate(unique_chroms)}
    
    # Prepare data for plotting
    y_positions = []
    x_values = []
    for chrom in unique_chroms:
        coords = chrom_coords[chrom]
        y_positions.extend([chrom_to_y[chrom]] * len(coords))
        x_values.extend(coords)
    
    # Determine the maximum chromosome length for scaling
    max_chrom_length = max([chrom_lengths[chrom] for chrom in unique_chroms])
    
    # Create the plot
    plt.figure(figsize=(12, max(6, len(unique_chroms) * 0.3)))
    plt.scatter(x_values, y_positions, color='blue', s=10, alpha=0.6)
    
    # Set y-axis labels to chromosome IDs
    plt.yticks(range(len(unique_chroms)), unique_chroms)
    
    # Adjust x-axis limits to the maximum chromosome length
    plt.xlim(0, max_chrom_length)
    
    # Draw vertical lines to indicate chromosome lengths
    for idx, chrom in enumerate(unique_chroms):
        chrom_length = chrom_lengths[chrom]
        plt.plot([chrom_length, chrom_length], [idx - 0.4, idx + 0.4], color='red', linestyle='--', linewidth=1)
    
    # Customize the plot
    plt.title('Scatter Plot of SNP Coordinates')
    plt.xlabel('Coordinate')
    plt.ylabel('Chromosome')
    plt.tight_layout()
    
    # Save the plot as a PDF
    plt.savefig(output_file, format='pdf')
    plt.close()
    print(f"Plot saved as {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python plot_snps.py <parental_snp_valid.vcf> <genome.fasta.fai>")
        sys.exit(1)
    vcf_file = sys.argv[1]
    fai_file = sys.argv[2]
    main(vcf_file, fai_file)
