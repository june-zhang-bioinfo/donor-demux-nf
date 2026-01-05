process DEMUXALOT {
    tag "${run_name}"
    publishDir "${params.demux_outdir}", mode: 'copy', pattern: "*.csv"

    input:
    tuple val(run_name), path(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam), path(barcodes), path(matrix_file)
    path demuxalot_script

    output:
    // Pass the matrix_file to the output tuple so Seurat can use it
    tuple val(run_name), path(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam), path("*.csv"), path(matrix_file), emit: demuxalot_out

    script:
    def samples_list_string = sample_names.join(',')
    """
    if [ ! -f "${tenx_bam}.bai" ]; then
        samtools index ${tenx_bam}
    fi
    python3 ${demuxalot_script} \
        --run_name "${run_name}" \
        --vcf_file "${filtered_vcf}" \
        --bam_file "${tenx_bam}" \
        --barcode_file "${barcodes}" \
        --sample_names "${samples_list_string}"
    """
}

process R_SEURAT_ANALYSIS {
    tag "${run_name}"

    input:
    // matrix_file is now staged directly in the work directory
    tuple val(run_name), path(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam), path(demux_results), path(matrix_file)
    path seurat_script

    output:
    path "PCA_plots.pdf", emit: pca_pdf
    path "SexGene_plots.pdf", emit: sex_gene_pdf

    script:
    def demux_csv = demux_results.find { it.name.endsWith('_likelihoods.csv') }
    def sex_list = sexes.collect{ "'$it'" }.join(',')
    def sample_list = sample_names.collect{ "'$it'" }.join(',')
    
    // Determine flag based on actual staged file extension
    def matrix_flag = matrix_file.name.endsWith('.h5') ? "--matrix_h5" : "--rds"
    
    """
    Rscript ${seurat_script} \
        --run_name "${run_name}" \
        --demux_csv "${demux_csv}" \
        ${matrix_flag} "${matrix_file}" \
        --sample_names "${sample_list}" \
        --sexes "${sex_list}"
    """
}

workflow SUBWF_DEMUX_ANALYSIS {
    take:
        demuxalot_input
        demuxalot_script
        seurat_script
    
    main:
        DEMUXALOT(
            demuxalot_input,
            demuxalot_script
        )
        
        R_SEURAT_ANALYSIS(
            DEMUXALOT.out.demuxalot_out,
            seurat_script
        )
    
}