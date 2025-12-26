#!/usr/bin/env python

import argparse
import os
import pandas as pd
from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps

def main():
    parser = argparse.ArgumentParser(description="Run Demuxalot on a 10x run using participant variants.")
    parser.add_argument("--run_name", required=True, help="The Nextflow run name (e.g., A1).")
    parser.add_argument("--vcf_file", required=True, help="Path to the filtered VCF file containing participant genotypes.")
    parser.add_argument("--bam_file", required=True, help="Path to the 10x BAM file for this run.")
    parser.add_argument("--sample_names", required=True, help="Comma-separated list of sample names (SM IDs).")
    
    args = parser.parse_args()

    # --- 1. Preparation ---
    sample_names = args.sample_names.strip("'").split(',')
    bam_path = args.bam_file
    barcode_path = os.path.join(os.path.dirname(bam_path), "barcodes.tsv")

    print(f"Run: {args.run_name}")
    print(f"Sample names: {sample_names}")
    print(f"BAM file: {bam_path}")
    print(f"Barcode file: {barcode_path}")

    # --- 2. Validate inputs ---
    if not os.path.exists(bam_path):
        print(f"ERROR: BAM file not found: {bam_path}")
        exit(1)
    if not os.path.exists(barcode_path):
        print(f"ERROR: Barcode file not found: {barcode_path}")
        exit(1)

    # --- 3. Load genotypes ---
    genotypes = ProbabilisticGenotypes(genotype_names=sample_names)
    genotypes.add_vcf(args.vcf_file)

    # --- 4. Load barcodes ---
    barcode_handler = BarcodeHandler.from_file(barcode_path)

    # --- 5. Count SNPs ---
    snps = count_snps(
        bamfile_location=bam_path,
        chromosome2positions=genotypes.get_chromosome2positions(),
        barcode_handler=barcode_handler
    )

    # --- 6. Predict posteriors ---
    likelihoods, posteriors = Demultiplexer.predict_posteriors(
        snps,
        genotypes=genotypes,
        only_singlets=True,
        barcode_handler=barcode_handler
    )

    # --- 7. Save results ---
    base_name = args.run_name 
    posteriors.to_csv(f'{base_name}_posteriors.csv')
    likelihoods.to_csv(f'{base_name}_likelihoods.csv')

    print(f"Demultiplexing complete for {args.run_name}!")

if __name__ == "__main__":
    main()