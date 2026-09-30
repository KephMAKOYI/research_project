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
process HAPLOTYPECALLER {
    scratch true

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id),
        path("${id}_hc.g.vcf.gz"),
        path("${id}_hc.g.vcf.gz.tbi"),
        emit: hc_gvcf
          
    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus} \
        -O ${id}_hc.g.vcf.gz
    """
}

//
workflow {
    // 1. Setup Input Channels
    ch_samples = Channel.fromPath(params.input).splitCsv(header: true).map { row -> tuple(row.sample, file(row.bam), file(row.bai)) }
    
    // 2. Linear Variant Calling
    hc_varcall = HAPLOTYPECALLER(ch_samples)
}
