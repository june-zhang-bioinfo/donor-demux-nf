#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// --- 1. Import Modules ---
if( !params.skip_vcf_prep ) {
    include { SUBWF_VCF_PREP } from './modules/vcf_preparation.nf'
}
include { SUBWF_DEMUX_ANALYSIS } from './modules/demultiplexing_analysis.nf'



// --- 2. Parameters ---
params.runs_csv = "/data/run_info.csv"
params.reference = "/data/GRCh38.primary_assembly.genome.fa"
params.tenx_dir = "/data/10x"
params.demux_outdir = "results/demuxalot"

// Conditional Demo Parameters (NEW)
params.skip_vcf_prep = false
params.vcf_path = "/data/text_dataset.vcf"
params.donor_ids = "113_113,349_350,352_353,39_39,40_40,41_41,42_42,43_43,465_466,596_597,597_598,632_633,633_634,660_661" // Comma-separated sample names for demo VCF
params.sexes = "Male,Male,Female,Female,Female,Male,Female,Female,Male,Female,Male,Male,Female,Female" // Comma-separated sexes for demo VCF

// Add this just after the parameter block
if (!params.skip_vcf_prep && !file(params.runs_csv).exists()) {
    exit 1, "Error: Full pipeline mode requires a valid run info CSV file at: ${params.runs_csv}"
}

// --- 3. Input Channels ---

// Parse CSV and group by run_name 
def runs_input_source = 
    (params.skip_vcf_prep || !file(params.runs_csv).exists()) ? 
    Channel.empty() : 
    Channel.fromPath(params.runs_csv)

// Apply the operators to the source, which will be empty in vcf skipping mode
runs_input_source
    .splitCsv(header:true)
    .map { row -> 
        tuple(row.run_name, row.bam_file, row.ID, row.LB, row.PL, row.PU, row.SM, row.sex)
    }
    .groupTuple()
    .set { runs_channel }

// Create channel for 10x BAMs
Channel.fromPath("${params.tenx_dir}/*.bam")
    .map { bam -> 
        def run_name = bam.parent.name
        def barcodes = file("${bam.parent}/barcodes.tsv")
        
        // Discover EITHER an .h5 or an .rds file in the same directory
        def matrix_files = files("${bam.parent}/*.{h5,rds}")
        def matrix = matrix_files.size() > 0 ? matrix_files[0] : file("NO_MATRIX_FILE")
        
        if( !barcodes.exists() ) {
            log.warn "Warning: Could not find barcodes.tsv for ${run_name}"
        }
        if( !matrix.exists() ) {
            log.warn "Warning: No .h5 or .rds found for ${run_name}"
        }
        
        return tuple(run_name, bam, barcodes, matrix) 
    }
    .set { tenx_bams_channel }


// --- 4. Conditional Workflow Execution ---
workflow {
    def final_vcf_channel
    
    if (!params.skip_vcf_prep) {
        vcf_prep_out = SUBWF_VCF_PREP(
            runs_channel,
            file(params.reference),
            file('scripts/change_read_group.py')
        )
        final_vcf_channel = vcf_prep_out.filtered_vcf
    } else {
        final_vcf_channel = tenx_bams_channel.map { run_name, bam, barcodes, matrix ->
            return tuple(
                run_name, 
                file(params.vcf_path), 
                params.donor_ids.split(','),
                params.sexes.split(',')
            )
        }
    }
    
    joined_channel = final_vcf_channel.join(tenx_bams_channel, by: 0)
    
    SUBWF_DEMUX_ANALYSIS(
        joined_channel,
        file('scripts/run_demuxalot.py'),
        file('scripts/seurat_analysis.R')
    )
}