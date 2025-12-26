#!/usr/bin/env python

import argparse
import pysam
import os

def main():
    parser = argparse.ArgumentParser(description="Change Read Group (RG) tags in BAM files using Pysam.")
    parser.add_argument("--run_name", required=True, help="The Nextflow run name.")
    parser.add_argument("--bam_files", required=True, help="Comma-separated list of input BAM file paths.")
    parser.add_argument("--ids", required=True, help="Comma-separated list of RG IDs.")
    parser.add_argument("--lbs", required=True, help="Comma-separated list of RG Library (LB) tags.")
    parser.add_argument("--pls", required=True, help="Comma-separated list of RG Platform (PL) tags.")
    parser.add_argument("--pus", required=True, help="Comma-separated list of RG Platform Unit (PU) tags.")
    parser.add_argument("--sms", required=True, help="Comma-separated list of RG Sample (SM) tags.")
    
    args = parser.parse_args()

    # --- 1. Parse Arguments ---
    
    # Split the comma-separated strings back into lists
    bam_files = args.bam_files.split(',')
    ids = args.ids.split(',')
    lbs = args.lbs.split(',')
    pls = args.pls.split(',')
    pus = args.pus.split(',')
    sms = args.sms.split(',')
    
    run_name = args.run_name
    output_dir = os.path.join('renamed_bams', run_name)
    os.makedirs(output_dir, exist_ok=True)

    # --- 2. Process BAM files ---
    
    for i, input_bam in enumerate(bam_files):
        output_bam = os.path.join(output_dir, os.path.basename(input_bam))
        
        # Define the new Read Group dictionary
        read_group = {
            'ID': ids[i],
            'LB': lbs[i],
            'PL': pls[i],
            'PU': pus[i],
            'SM': sms[i]
        }
        
        try:
            with pysam.AlignmentFile(input_bam, "rb") as infile:
                # Update the header dictionary
                header = infile.header.to_dict()
                # Overwrite or create the 'RG' tag list with the new group
                header['RG'] = [read_group]
                
                with pysam.AlignmentFile(output_bam, "wb", header=header) as outfile:
                    for read in infile.fetch():
                        # Set the RG tag for the individual read
                        read.set_tag('RG', read_group['ID'])
                        outfile.write(read)
            
            print(f"Processed: {output_bam}")
            
        except Exception as e:
            print(f"ERROR processing {input_bam}: {e}")
            exit(1)

if __name__ == "__main__":
    main()