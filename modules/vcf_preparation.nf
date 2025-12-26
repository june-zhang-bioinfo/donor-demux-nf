// modules/vcf_preparation.nf

// Step 1: Change Read Group for participant BAMs
process CHANGE_READ_GROUP {
    container '/PHShome/zz005/docker/donor_demux.sif'
    
    input:
    tuple val(run_name), val(bam_files), val(ids), val(lbs), val(pls), val(pus), val(sms), val(sexes)
    path change_read_group_script
    
    output:
    tuple val(run_name), path("renamed_bams/${run_name}/*.bam"), val(sms), val(sexes), emit: renamed_bams
    
    script:
    def bam_string = bam_files.join(',')
    def id_string = ids.join(',')
    def lb_string = lbs.join(',')
    def pl_string = pls.join(',')
    def pu_string = pus.join(',')
    def sm_string = sms.join(',')
    
    """
    python3 ${change_read_group_script} \\
        --run_name "${run_name}" \\
        --bam_files "${bam_string}" \\
        --ids "${id_string}" \\
        --lbs "${lb_string}" \\
        --pls "${pl_string}" \\
        --pus "${pu_string}" \\
        --sms "${sm_string}"
    """
}

// Step 2: FreeBayes variant calling on participant BAMs
process FREEBAYES_CALL {
    container '/PHShome/zz005/docker/donor_demux.sif'
    
    input:
    tuple val(run_name), path(bam_files), val(sms), val(sexes)
    path reference_genome
    
    output:
    tuple val(run_name), path("freebayes/${run_name}/filtered.vcf"), val(sms), val(sexes), emit: filtered_vcf
    
    script:
    """
    mkdir -p freebayes/${run_name}
    
    # Create BAM list file
    bam_list="freebayes/${run_name}/bam_list.txt"
    for f in ${bam_files}; do 
        echo \$(readlink -f \$f) >> \$bam_list
    done
    
    # Index BAM files if not already indexed
    for f in ${bam_files}; do
        if [ ! -f "\${f}.bai" ];
            samtools index \$f
        fi
    done
    
    # Run FreeBayes
    freebayes -f ${reference_genome} --bam-list \$bam_list > freebayes/${run_name}/variants.vcf
    
    # Filter variants
    bcftools view -i 'QUAL>30 & INFO/DP>10 & INFO/AF>0.1' \\
        freebayes/${run_name}/variants.vcf > freebayes/${run_name}/filtered.vcf
    """
}

// Subworkflow to chain these two processes
workflow SUBWF_VCF_PREP {
    take: 
    runs_channel
    reference_genome 
    read_group_script
    
    emit: filtered_vcf
    
    CHANGE_READ_GROUP(runs_channel, read_group_script) // <-- Pass the script here
    
    FREEBAYES_CALL(
        CHANGE_READ_GROUP.out.renamed_bams,
        reference_genome // <-- Pass the reference here
    )
    
    emit: FREEBAYES_CALL.out.filtered_vcf
}