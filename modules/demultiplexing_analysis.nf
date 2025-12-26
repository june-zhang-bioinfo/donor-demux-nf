process DEMUXALOT {
    container '/PHShome/zz005/docker/donor_demux.sif'
    publishDir "${params.demux_outdir}/${run_name}", mode: 'copy'

    input:
    tuple val(run_name), path(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam)
    path demuxalot_script

    output:
    tuple val(run_name), val(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam), path("*.csv"), emit: demuxalot_out

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
        --sample_names "${samples_list_string}"
    """
}

process R_SEURAT_ANALYSIS {
    container 'r_eval:1.0'
    publishDir "${params.demux_outdir}/${run_name}", mode: 'copy'

    input:
    tuple val(run_name), val(filtered_vcf), val(sample_names), val(sexes), path(tenx_bam), path(demux_results)
    path seurat_script

    output:
    path "${run_name}_PCA_plots.pdf", emit: pca_pdf
    path "${run_name}_SexGene_plots.pdf", emit: sex_gene_pdf

    script:
    def matrix_dir = tenx_bam.parent.parent
    def rds_files = file("${matrix_dir}/*.rds")
    def rds_file = rds_files.size() > 0 ? rds_files[0] : null
    
    def h5_path = "${matrix_dir}/filtered_feature_bc_matrix.h5"
    def h5_file = (!rds_file && file(h5_path).exists()) ? file(h5_path) : null
    
    def demux_csv = demux_results.find { it.name.endsWith('_posteriors.csv') }
    def sex_list = sexes.collect{ "'$it'" }.join(',')
    def sample_list = sample_names.collect{ "'$it'" }.join(',')
    
    def r_cmd = ["Rscript ${seurat_script}",
                 "--run_name '${run_name}'",
                 "--demux_csv '${demux_csv}'",
                 "--sample_names '${sample_list}'",
                 "--sexes '${sex_list}'"]
    
    if (rds_file) {
        r_cmd << "--rds '${rds_file}'"
    } else if (h5_file) {
        r_cmd << "--matrix_h5 '${h5_file}'"
    } else {
        error "No RDS or H5 file found in ${matrix_dir}"
    }
    
    """
    ${r_cmd.join(' ')}
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
    
    // NO EMIT BLOCK HERE - Remove it completely!
}