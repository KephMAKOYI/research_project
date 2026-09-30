#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process showJavaVersion {
    publishDir "${params.outdir}/", mode: 'copy', overwrite: false
    output:
    file 'java_version.txt'

    """
    java -version 2> java_version.txt
    """
}


//For DeepVariant
process GLNEXUS_MERGE {
    //publishDir "${params.output}/vcf", mode: 'copy', overwrite: true    

    input:
    path(gvcfs_and_indices) // A list of all .g.vcf.gz files

    output:
    tuple path("dv_cohort.vcf.gz"),
          path("dv_cohort.vcf.gz.tbi"),
          emit: join_gvcf

    script:
    // Filter the list to find only the .vcf.gz files (indices are just ignored but present)
    def variants = gvcfs_and_indices
        .findAll { it.name.endsWith(".vcf.gz") }
        .collect { "${it}" }
        .join(" ")
		
    """
    # Use the DeepVariantWGS preset for best results with your 26GB files
    glnexus_cli \
        --config DeepVariantWGS \
        --threads ${task.cpus} \
        ${variants} | bcftools view -Oz -o dv_cohort.vcf.gz

    tabix -p vcf dv_cohort.vcf.gz
    """
}

//
process SELECTVARIANTS {
    //publishDir "${params.output}/", mode: 'copy', overwrite: true

    input:
    tuple path(joint_vcf), path(joint_tbi)

    output:
    tuple path("dv_raw_snps.vcf.gz"),
          path("dv_raw_snps.vcf.gz.tbi"),
          emit: snps
    tuple path("dv_raw_indels.vcf.gz"),
          path("dv_raw_indels.vcf.gz.tbi"),
          emit: indels

    script:
    """
    #
    gatk SelectVariants \
        -R ${params.reference_files.ref} \
        -V ${joint_vcf} \
        --select-type SNP \
        -O dv_raw_snps.vcf.gz

    gatk SelectVariants \
        -R ${params.reference_files.ref} \
        -V ${joint_vcf} \
        --select-type INDEL \
        -O dv_raw_indels.vcf.gz
    """
}

//
process VARIANTFILTRATION {
    //publishDir "${params.output}/", mode: 'copy', overwrite: true
    input:
    tuple path(snps), path(snps_tbi)
    tuple path(indels), path(indels_tbi)

    output:
    tuple path("dv_raw_filtered_snps.vcf.gz"),
          path("dv_raw_filtered_snps.vcf.gz.tbi"),
          emit: filter_snps
    tuple path("dv_raw_filtered_indels.vcf.gz"),
          path("dv_raw_filtered_indels.vcf.gz.tbi"),
          emit: filter_indels

    script:
    """
    # Filter SNPs
    bcftools filter -S . -e 'FMT/DP<10 || FMT/GQ<20' ${snps} | \
    bcftools view -i 'QUAL>30 && FILTER="PASS"' -O z -o dv_raw_filtered_snps.vcf.gz
    bcftools index -t dv_raw_filtered_snps.vcf.gz

    # Filter INDELS
    bcftools filter -S . -e 'FMT/DP<10 || FMT/GQ<20' ${indels} | \
    bcftools view -i 'QUAL>50 && FILTER="PASS"' -O z -o dv_raw_filtered_indels.vcf.gz
    bcftools index -t dv_raw_filtered_indels.vcf.gz
    """
}

//
process MERGE_VCF {
    //publishDir "${params.outdir}/vcf", mode: 'copy'
    input:
    tuple path(snp), path(stbi)
    tuple path(indel), path(itbi)

    output:
    tuple path("dv_joint_filtered.vcf.gz"),
          path("dv_joint_filtered.vcf.gz.tbi"),
          emit: vcf_filtered_concat

    script:
    """
    # Concatenate and index
    bcftools concat -a ${snp} ${indel} -O z -o dv_joint_filtered.vcf.gz
    bcftools index -t dv_joint_filtered.vcf.gz
    """
}

//
process VCF_QC {
    //publishDir "${params.outdir}/vcf_qc", mode: 'copy'
    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("dv_joint_sorted.vcf.gz"),
          path("dv_joint_sorted.vcf.gz.tbi")
    tuple path("dv_joint_normalized.vcf.gz"), 
          path("dv_joint_normalized.vcf.gz.tbi")
    tuple path("dv_joint_normalized_stats.txt")

    script:
    """
    # 1. Sort the VCF (if not already sorted)
    bcftools sort ${vcf} -O z -o dv_joint_sorted.vcf.gz
    bcftools index -t dv_joint_sorted.vcf.gz

    # 2. Decompose multi-allelics (-m -), normalize, and filter for Left-Alignment/Trim (-f)
    bcftools norm -m -both -f ${params.reference_files.ref} dv_joint_sorted.vcf.gz -O z -o dv_joint_normalized.vcf.gz
    bcftools index -t dv_joint_normalized.vcf.gz

    # Run bcftools stats
    bcftools stats dv_joint_normalized.vcf.gz > dv_joint_normalized_stats.txt
    """
}

process VEP {
    scratch true   // ⭐ensures node-local temp cleanup
    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("dv_vep_annotated_final.vcf.gz"),
          path("dv_vep_annotated_final.vcf.gz.tbi"),
          path("dv_vep_annotated_vep_report.html")

    script:
    """
    vep -i ${vcf} \
        -o dv_vep_annotated_final.vcf.gz \
        --vcf \
        --stats_file dv_vep_annotated_vep_report.html \
        --compress_output bgzip \
        --assembly ${params.assembly} \
        --cache \
        --dir_cache ${params.vep_cache_dir} \
        --plugin AlphaMissense,file=${params.alphamissense},cols=all \
        --force_overwrite \
        --species homo_sapiens \
        --offline \
        --everything \
        --fork ${task.cpus}

    tabix -p vcf dv_vep_annotated_final.vcf.gz
    """
}

// WORKFLOW

workflow {
    // This finds pairs of files (vcf.gz AND vcf.gz.tbi) and stages them together
    
    // -----------------------------------------------------------
    // STEP 1. Channel setup (same as yours)
    // -----------------------------------------------------------

    //gvcf_ch = Channel.fromPath(params.dv_gvcfs).map { file -> tuple(file, file + ".tbi") }.collect()

    gvcf_ch = Channel.fromPath(params.input).splitCsv(header: true).map { row -> tuple(file(row.gvcf), file(row.tbi)) }.collect()

    // -----------------------------------------------------------
    // STEP 2. Combine all GVCFs into one cohort GVCF
    // -----------------------------------------------------------
    
    glnexus_ch = GLNEXUS_MERGE(gvcf_ch).join_gvcf
    ch_select_variant = SELECTVARIANTS(glnexus_ch)
    
    // -----------------------------------------------------------
    // STEP 3 — VariantFiltration
    // -----------------------------------------------------------
    
    filter_ch = VARIANTFILTRATION(ch_select_variant.snps,ch_select_variant.indels)
    ch_concat_vcf = MERGE_VCF(filter_ch.filter_snps, filter_ch.filter_indels)
    VCF_QC(ch_concat_vcf)
    
    // -----------------------------------------------------------
    // STEP 4 — VEP annotation
    // -----------------------------------------------------------
    VEP(ch_concat_vcf)
}
