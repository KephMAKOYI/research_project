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

//
process CombineGVCFs {
    input:
    path gvcfs_and_indices

    output:
    tuple path("cohort.g.vcf.gz"),
          path("cohort.g.vcf.gz.tbi"),
          emit: combined_gvcf
    
    script:

    def variants = gvcfs_and_indices.findAll { it.name.endsWith(".vcf.gz") }.collect { "-V ${it}" }.join(" ")

    """
    gatk CombineGVCFs \
        -R ${params.reference_files.ref} \
        ${variants} \
        -O cohort.g.vcf.gz
    """
}

//
process GenotypeGVCFs {
    input:
    tuple path(cohort_vcf), path(cohort_tbi)

    output:
    tuple path("hc_genotype.vcf.gz"),
          path("hc_genotype.vcf.gz.tbi"),
          emit: vcf_shard

    script:
    """
    gatk GenotypeGVCFs \
        -R ${params.reference_files.ref} \
        -V ${cohort_vcf} \
        -O hc_genotype.vcf.gz
    """
}

// ------------------------------------------------------------------------
// VQSR: SNP Mode
// ------------------------------------------------------------------------
process VQSR_SNP {
    input:
    tuple path(genotype), path(tbis)

    output:
    path("hc_genotype_vqsr_snp.recal"), emit: recal
    path("hc_genotype_vqsr_snp.recal.idx"), emit: recal_idx
    path("hc_genotype_vqsr_snp.tranches"), emit: tranches

    script:
    """
    gatk VariantRecalibrator \
        -R ${params.reference_files.ref} \
        -V ${genotype} \
        -AS \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.reference_files.hapmap} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.reference_files.omni} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.reference_files.phase1} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.reference_files.dbsnp} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O hc_genotype_vqsr_snp.recal \
        --tranches-file hc_genotype_vqsr_snp.tranches

    gatk IndexFeatureFile -I hc_genotype_vqsr_snp.recal
    """
}

// ------------------------------------------------------------------------
// VQSR: INDEL Mode
// ------------------------------------------------------------------------
process VQSR_INDEL {
    input:
    tuple path(genotype), path(tbis)

    output:
    path("hc_genotype_vqsr_indel.recal"), emit: recal
    path("hc_genotype_vqsr_indel.recal.idx"), emit: recal_idx
    path("hc_genotype_vqsr_indel.tranches"), emit: tranches

    script:
    """
    gatk VariantRecalibrator \
        -R ${params.reference_files.ref} \
        -V ${genotype} \
        -AS \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${params.reference_files.mills} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.reference_files.dbsnp} \
        --resource:known_indels,known=false,training=true,truth=false,prior=10.0 ${params.reference_files.known_indels} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode INDEL \
        -O hc_genotype_vqsr_indel.recal \
        --tranches-file hc_genotype_vqsr_indel.tranches

    gatk IndexFeatureFile -I hc_genotype_vqsr_indel.recal
    """
}

// ------------------------------------------------------------------------
// Apply VQSR (Apply both SNP and INDEL filters)
// ------------------------------------------------------------------------
process ApplyVQSR {
    scratch true   // ⭐ensures node-local temp cleanup        

    input:
    tuple path(genotype), path(tbi)
    path(snp_recal)
    path(snp_idx)
    path(snp_tranches)
    path(indel_recal)
    path(indel_idx)
    path(indel_tranches)

    output:
    tuple path("hc_final_recalibrated.vcf.gz"),
          path("hc_final_recalibrated.vcf.gz.tbi"),
          emit: vcf

    script:
    """
    set -euo pipefail

    # 1. Apply SNP Recalibration
    gatk ApplyVQSR \
        -R ${params.reference_files.ref} \
        -V ${genotype} \
        -mode SNP \
        --recal-file ${snp_recal} \
        --tranches-file ${snp_tranches} \
        --truth-sensitivity-filter-level 99.0 \
        -O hc_tmp_snp_recal.vcf.gz

    # 2. Apply INDEL Recalibration (takes output of step 1 as input)
    gatk ApplyVQSR \
        -R ${params.reference_files.ref} \
        -V hc_tmp_snp_recal.vcf.gz \
        -mode INDEL \
        --recal-file ${indel_recal} \
        --tranches-file ${indel_tranches} \
        --truth-sensitivity-filter-level 99.0 \
        -O hc_final_recalibrated.vcf.gz
    """
}

//
process SELECTVARIANT {
    publishDir "${params.outdir}/${params.hc_joint}", mode: 'copy', overwrite: false

    input:
    tuple path(variant_recal), path(var_tbi)
  
    output:
    tuple path("hc_joint_variant.vcf.gz"),
          path("hc_joint_variant.vcf.gz.tbi"),
          emit: var_select

    script:
    """
    gatk SelectVariants \
        --exclude-filtered \
        -V ${variant_recal} \
        -O hc_joint_variant.vcf.gz
    """
}

//
process SELECTVARIANTS {
    //publishDir "${params.output}/", mode: 'copy', overwrite: true

    input:
    tuple path(recal_vcf), path(recal_tbi)

    output:
    tuple path("hc_raw_snps.vcf.gz"),
          path("hc_raw_snps.vcf.gz.tbi"),
          emit: snps
    tuple path("hc_raw_indels.vcf.gz"),
          path("hc_raw_indels.vcf.gz.tbi"),
          emit: indels

    script:
    """
    #
    gatk SelectVariants \
        -R ${params.reference_files.ref} \
        -V ${recal_vcf} \
        --select-type SNP \
        -O hc_raw_snps.vcf.gz

    gatk SelectVariants \
        -R ${params.reference_files.ref} \
        -V ${recal_vcf} \
        --select-type INDEL \
        -O hc_raw_indels.vcf.gz
    """
}

//
process GATK_FILTRATION {
    //publishDir "${params.output}/", mode: 'copy', overwrite: true

    input:
    tuple path(snps), path(snps_tbi)
    tuple path(indels), path(indels_tbi)

    output:
    tuple path("gatk_hc_raw_filtered_snps.vcf.gz"),
          path("gatk_hc_raw_filtered_snps.vcf.gz.tbi"),
          emit: gatk_filter_snps
    tuple path("gatk_hc_raw_filtered_indels.vcf.gz"),
          path("gatk_hc_raw_filtered_indels.vcf.gz.tbi"),
          emit: gatk_filter_indels

    script:
    """
    # Filter SNPs
    gatk VariantFiltration \
        -R ${params.reference_files.ref} \
        -V ${snps} \
        -O gatk_hc_raw_filtered_snps.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
        -genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 10"

    # Filter INDELS
    gatk VariantFiltration \
        -R ${params.reference_files.ref} \
        -V ${indels} \
        -O gatk_hc_raw_filtered_indels.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10" \
        -genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 10"
    """
}

//
process BCF_FILTRATION {
    //publishDir "${params.output}/", mode: 'copy', overwrite: true

    input:
    tuple path(snps), path(snps_tbi)
    tuple path(indels), path(indels_tbi)

    output:
    tuple path("bcf_hc_raw_filtered_snps.vcf.gz"),
          path("bcf_hc_raw_filtered_snps.vcf.gz.tbi"),
          emit: bcf_filter_snps
    tuple path("bcf_hc_raw_filtered_indels.vcf.gz"),
          path("bcf_hc_raw_filtered_indels.vcf.gz.tbi"),
          emit: bcf_filter_indels

    script:
    """
    # Filter SNPs
    bcftools filter -S . -e 'FMT/DP<10 || FMT/GQ<20' ${snps} | \
    bcftools view -i 'AC>=1 && F_MISSING<=0.2 && QUAL>30 && AF>=0.05' -O z -o bcf_hc_raw_filtered_snps.vcf.gz
    bcftools index -t bcf_hc_raw_filtered_snps.vcf.gz

    # Filter INDELS
    bcftools filter -S . -e 'FMT/DP<10 || FMT/GQ<20' ${indels} | \
    bcftools view -i 'AC>=1 && F_MISSING<=0.2 && QUAL>50 && AF>=0.05' -O z -o bcf_hc_raw_filtered_indels.vcf.gz
    bcftools index -t bcf_hc_raw_filtered_indels.vcf.gz
    """
}

//
process GATK_MERGE_VCF {
    //publishDir "${params.outdir}/${params.hc_joint}", mode: 'copy', overwrite: true

    input:
    tuple path(snp), path(stbi)
    tuple path(indel), path(itbi)

    output:
    tuple path("gatk_hc_joint_filtered.vcf.gz"),
          path("gatk_hc_joint_filtered.vcf.gz.tbi"),
          emit: gatk_vcf_filtered_concat

    script:
    """
    # Concatenate and index
    bcftools concat -a ${snp} ${indel} -O z -o gatk_hc_joint_filtered.vcf.gz
    bcftools index -t gatk_hc_joint_filtered.vcf.gz
    """
}

//
process BCF_MERGE_VCF {
    //publishDir "${params.outdir}/${params.hc_joint}", mode: 'copy', overwrite: true

    input:
    tuple path(snp), path(stbi)
    tuple path(indel), path(itbi)

    output:
    tuple path("bcf_hc_joint_filtered.vcf.gz"),
          path("bcf_hc_joint_filtered.vcf.gz.tbi"),
          emit: bcf_vcf_filtered_concat

    script:
    """
    # Concatenate and index
    bcftools concat -a ${snp} ${indel} -O z -o bcf_hc_joint_filtered.vcf.gz
    bcftools index -t bcf_hc_joint_filtered.vcf.gz
    """
}

//
process GATK_VCF_QC {
    //publishDir "${params.outdir}/${params.vcf_qc}", mode: 'copy', overwrite: true

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("gatk_hc_joint_sorted.vcf.gz"),
          path("gatk_hc_joint_sorted.vcf.gz.tbi")
    tuple path("gatk_hc_joint_normalized.vcf.gz"),
          path("gatk_hc_joint_normalized.vcf.gz.tbi"),
          emit: gatk_joint_norm
    tuple path("gatk_hc_joint_normalized_stats.txt")

    script:
    """
    # 1. Sort the VCF (if not already sorted)
    bcftools sort ${vcf} -O z -o gatk_hc_joint_sorted.vcf.gz
    bcftools index -t gatk_hc_joint_sorted.vcf.gz

    # 2. Decompose multi-allelics (-m -), normalize, and filter for Left-Alignment/Trim (-f)
    bcftools norm -m -both -f ${params.reference_files.ref} gatk_hc_joint_sorted.vcf.gz -O z -o gatk_hc_joint_normalized.vcf.gz
    bcftools index -t gatk_hc_joint_normalized.vcf.gz

    # Run bcftools stats
    bcftools stats gatk_hc_joint_normalized.vcf.gz > gatk_hc_joint_normalized_stats.txt
    """
}

//
process BCF_VCF_QC {
    //publishDir "${params.outdir}/${params.vcf_qc}", mode: 'copy', overwrite: true

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("bcf_hc_joint_sorted.vcf.gz"),
          path("bcf_hc_joint_sorted.vcf.gz.tbi")
    tuple path("bcf_hc_joint_normalized.vcf.gz"),
          path("bcf_hc_joint_normalized.vcf.gz.tbi"),
          emit : bcf_joint_norm
    tuple path("bcf_hc_joint_normalized_stats.txt")

    script:
    """
    # 1. Sort the VCF (if not already sorted)
    bcftools sort ${vcf} -O z -o bcf_hc_joint_sorted.vcf.gz
    bcftools index -t bcf_hc_joint_sorted.vcf.gz

    # 2. Decompose multi-allelics (-m -), normalize, and filter for Left-Alignment/Trim (-f)
    bcftools norm -m -both -f ${params.reference_files.ref} bcf_hc_joint_sorted.vcf.gz -O z -o bcf_hc_joint_normalized.vcf.gz
    bcftools index -t bcf_hc_joint_normalized.vcf.gz

    # Run bcftools stats
    bcftools stats bcf_hc_joint_normalized.vcf.gz > bcf_hc_joint_normalized_stats.txt
    """
}


//
process GATK_VEP {
    scratch true   // ⭐ensures node-local temp cleanup
    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("gatk_vep_annotated_final.vcf.gz"),
          path("gatk_vep_annotated_final.vcf.gz.tbi"),
          path("gatk_vep_annotated_vep_report.html")

    script:
    """
    vep -i ${vcf} \
        -o gatk_vep_annotated_final.vcf.gz \
        --vcf \
        --stats_file gatk_vep_annotated_vep_report.html \
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

    tabix -p vcf gatk_vep_annotated_final.vcf.gz
    """
}


//
process BCF_VEP {
    scratch true   // ⭐ensures node-local temp cleanup
    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path("bcf_vep_annotated_final.vcf.gz"),
          path("bcf_vep_annotated_final.vcf.gz.tbi"),
          path("bcf_vep_annotated_vep_report.html")

    script:
    """
    vep -i ${vcf} \
        -o bcf_vep_annotated_final.vcf.gz \
        --vcf \
        --stats_file bcf_vep_annotated_vep_report.html \
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

    tabix -p vcf bcf_vep_annotated_final.vcf.gz
    """
}


// WORKFLOW

workflow {
    // ch_chroms = Channel.fromList((1..22).collect{"chr$it"} + ['chrX','chrY','chrM'])
    
    //gvcf_ch = Channel.fromPath(params.hc_gvcfs).map { file -> tuple(file, file + ".tbi") }.collect()
    gvcf_ch = Channel.fromPath(params.input).splitCsv(header: true).map { row -> tuple(file(row.vcf), file(row.tbi)) }.collect()
    
    // Split into separate lists
    combined_gvcf_ch = CombineGVCFs(gvcf_ch).combined_gvcf
    
    // Genotype per chromosome
    genotype_shards = GenotypeGVCFs(combined_gvcf_ch).vcf_shard
    
    // 2. RUN the recalibrator processes (This was the missing step)
    snp_recal_ch = VQSR_SNP(genotype_shards)
    indel_recal_ch = VQSR_INDEL(genotype_shards)
    
    // 3. Pass the outputs of those runs into ApplyVQSR
    final_vcf_ch = ApplyVQSR(
       genotype_shards,
       snp_recal_ch.recal,
       snp_recal_ch.recal_idx,
       snp_recal_ch.tranches,
       indel_recal_ch.recal,
       indel_recal_ch.recal_idx,
       indel_recal_ch.tranches
    ).vcf
    
    // 4. Select PASS
    select_var_ch = SELECTVARIANT(final_vcf_ch).var_select
    
    // 5.
    ch_select_variant = SELECTVARIANTS(final_vcf_ch)
    //
    ch_variant_filtration_gatk = GATK_FILTRATION(ch_select_variant.snps,ch_select_variant.indels)
    
    ch_concat_vcf_gatk = GATK_MERGE_VCF(ch_variant_filtration_gatk.gatk_filter_snps, ch_variant_filtration_gatk.gatk_filter_indels)
    
    ch_vcf_qc_gatk = GATK_VCF_QC(ch_concat_vcf_gatk).gatk_joint_norm
    
    GATK_VEP(ch_vcf_qc_gatk)
    
    ch_variant_filtration_bcf = BCF_FILTRATION(ch_select_variant.snps,ch_select_variant.indels)
    
    ch_concat_vcf_bcf = BCF_MERGE_VCF(ch_variant_filtration_bcf.bcf_filter_snps, ch_variant_filtration_bcf.bcf_filter_indels)
    
    ch_vcf_qc_bcf = BCF_VCF_QC(ch_concat_vcf_bcf).bcf_joint_norm
    
    BCF_VEP(ch_vcf_qc_bcf)
}
