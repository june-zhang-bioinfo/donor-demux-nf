// modules/vcf_preparation.nf

process CHANGE_READ_GROUP {
    tag "${run_name}"

    input:
    tuple val(run_name), path(bam_files), val(donor_ids), val(lbs), val(pls), val(pus), val(sms), val(sexes)
    path change_read_group_script
    
    output:
    // Bundling .bam and .bai ensures they are staged together for FreeBayes
    tuple val(run_name), path("renamed_bams/${run_name}/*.{bam,bai}"), val(donor_ids), val(sexes), emit: renamed_bams
    
    script:
    def bam_string = bam_files.join(',')
    def id_string = donor_ids.join(',')
    def lb_string = lbs.join(',')
    def pl_string = pls.join(',')
    def pu_string = pus.join(',')
    def sm_string = sms.join(',')
    
    """
    # 1. Index input BAMs for pysam.fetch()
    for f in *.bam; do
        if [ ! -f "\${f}.bai" ]; then samtools index \$f; fi
    done

    # 2. Run the renaming script (using your robust .tmp + rename version)
    mkdir -p renamed_bams/${run_name}
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

process FREEBAYES_CALL {
    tag "${run_name}"
    
    input:
    tuple val(run_name), path(bam_files), val(donor_ids), val(sexes)
    path reference_genome
    
    output:
    tuple val(run_name), path("freebayes/filtered.vcf"), val(donor_ids), val(sexes), emit: filtered_vcf
    
    script:
    """
    mkdir -p freebayes
    
    # Generate BAM list for FreeBayes
    ls *.bam > bam_list.txt
    
    # Run FreeBayes
    # Note: FreeBayes will merge samples that share the same SM tag in the header
    freebayes -f ${reference_genome} --bam-list bam_list.txt > freebayes/variants.vcf
    
    # Filter variants
    bcftools view -i 'QUAL>30 & INFO/DP>10 & INFO/AF>0.1' \\
        freebayes/variants.vcf > freebayes/filtered.vcf
    """
}

workflow SUBWF_VCF_PREP {
    take: 
        runs_channel
        reference_genome 
        read_group_script
    
    main:
        CHANGE_READ_GROUP(runs_channel, read_group_script)
        
        FREEBAYES_CALL(
            CHANGE_READ_GROUP.out.renamed_bams,
            reference_genome
        )
    
    emit:
        filtered_vcf = FREEBAYES_CALL.out.filtered_vcf.map { run_name, vcf, donor_ids, sexes ->
            // 1. Pair sample names with their sexes
            // 2. Remove duplicates
            // 3. Sort numerically/alphabetically by sample name (it[0])
            def paired = [donor_ids, sexes].transpose().unique().sort { it[0] }

            // Extract the unique, sorted lists back
            def unique_donor_ids = paired.collect { it[0] }
            def unique_sexes = paired.collect { it[1] }
            
            return tuple(run_name, vcf, unique_donor_ids, unique_sexes)
        }
}