#!/usr/bin/env python
import argparse
import pysam
import os
import sys

def change_read_group(input_bam, output_bam, read_group):
    """
    Change read groups in a BAM file and index it.
    """
    temp_bam = output_bam + ".tmp"
    
    try:
        # Open the input BAM file
        with pysam.AlignmentFile(input_bam, "rb") as infile:
            # Update the header with the new read group
            header = infile.header.to_dict()
            header['RG'] = [read_group]
            
            # Write to temporary file first
            with pysam.AlignmentFile(temp_bam, "wb", header=header) as outfile:
                for read in infile.fetch(until_eof=True):
                    # Change the read group
                    read.set_tag('RG', read_group['ID'])
                    outfile.write(read)
                os.sync()
        
        # Only rename if write was successful
        os.rename(temp_bam, output_bam)
        
        # Index the file
        pysam.index(output_bam)
        
        print(f"✓ Processed {input_bam} -> {output_bam}")
        
    except Exception as e:
        # Clean up temp file if something failed
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
        print(f"✗ ERROR processing {input_bam}: {e}", file=sys.stderr)
        raise

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_name", required=True)
    parser.add_argument("--bam_files", required=True)
    parser.add_argument("--ids", required=True)
    parser.add_argument("--lbs", required=True)
    parser.add_argument("--pls", required=True)
    parser.add_argument("--pus", required=True)
    parser.add_argument("--sms", required=True)
    args = parser.parse_args()

    # Split the comma-separated strings passed by Nextflow
    bam_list = args.bam_files.split(',')
    id_list = args.ids.split(',')
    lb_list = args.lbs.split(',')
    pl_list = args.pls.split(',')
    pu_list = args.pus.split(',')
    sm_list = args.sms.split(',')

    output_dir = os.path.join('renamed_bams', args.run_name)
    os.makedirs(output_dir, exist_ok=True)

    for i in range(len(bam_list)):
        input_bam = bam_list[i]
        output_bam = os.path.join(output_dir, os.path.basename(input_bam))
        
        read_group = {
            'ID': id_list[i],
            'LB': lb_list[i],
            'PL': pl_list[i],
            'PU': pu_list[i],
            'SM': sm_list[i]
        }
        
        change_read_group(input_bam, output_bam, read_group)

if __name__ == "__main__":
    main()