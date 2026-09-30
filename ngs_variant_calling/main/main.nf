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

// Define the input path (adjust the path to your actual directory)
process FASTQC {
    publishDir "${params.outdir}/${params.fastqc_report}", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(read_1), path(read_2)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc -t ${task.cpus} -q ${read_1} ${read_2}
    """
}

//
process TRIMGALORE {
    publishDir "${params.outdir}/trim", mode: 'copy', overwrite: false

    input:
    tuple val(id), val(lane), path(r1), path(r2)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_R1.trim.fq.gz"),
          path("${id}_${lane}_R2.trim.fq.gz"),
          emit: trimmed_fq

    script:
    """
    trim_galore \
        --paired \
        --illumina \
        -j ${task.cpus} \
        --basename ${id}_${lane} \
        $r1 $r2

    mv ${id}_${lane}_val_1.fq.gz ${id}_${lane}_R1.trim.fq.gz
    mv ${id}_${lane}_val_2.fq.gz ${id}_${lane}_R2.trim.fq.gz
    """
}

//
process BWA {
    scratch true   // ⭐ensures node-local temp cleanup

    input:
    tuple val(id), val(lane), path(r1), path(r2)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_markdup.bam"),
          path("${id}_${lane}_markdup.bam.bai"),
          emit: markdup_bam

    script:
    """
    set -o pipefail

    # 1. Alignment against reference genome
    bwa mem -t ${task.cpus} \
       ${params.reference_files.ref} \
       $r1 $r2 | \
    samtools sort -n -@ ${task.cpus} -m 4G -u -T ${params.tmp_dir}/${id}_${lane}_nsort - | \
    samtools fixmate -m -u -@ 8 - - | \
    samtools sort -@ ${task.cpus} -m 4G -T ${params.tmp_dir}/${id}_${lane}_csort - | \
    samtools markdup -@ ${task.cpus} - - | \
    samtools addreplacerg \
       -r '@RG\\tID:${id}.${lane}\\tSM:${id}\\tLB:lib1\\tPL:ILLUMINA\\tPU:${lane}' \
       -o ${id}_${lane}_markdup.bam -
	
    # 2. Indexing BAM file
    samtools index ${id}_${lane}_markdup.bam
    """
}

//
process BASERECALIBRATOR {
    scratch true   // ⭐ensures node-local temp cleanup

    input:
    tuple val(id), val(lane), path(bam), path(bai)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}.table"), 
          emit: table
    
    script:
    """
    gatk BaseRecalibrator \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        --known-sites ${params.reference_files.dbsnp} \
        -O ${id}_${lane}.table
    """
}

//
process ApplyBQSR {
    input:
    tuple val(id), val(lane), path(bam), path(bai), path(table)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_bqsr.bam"),
          path("${id}_${lane}_bqsr.bam.bai"),
          emit: bqsr_bam

    script:
    """
    gatk ApplyBQSR \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        --bqsr-recal-file ${table} \
        -O ${id}_${lane}_bqsr.bam

    samtools index ${id}_${lane}_bqsr.bam
    """
}


//
process MERGE_BAMS {
    publishDir "${params.outdir}/${params.bam}", mode: 'copy', overwrite: false
    scratch true

    input:
    tuple val(id), path(bams)

    output:
    tuple val(id),
          path("${id}_merged.bam"),
          path("${id}_merged.bam.bai"),
          emit: merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} \
          -o ${id}_merged.bam \
             ${bams.join(' ')}
    samtools index ${id}_merged.bam
    """
}

//
process BAM_QC {
    publishDir "${params.outdir}/${params.bam_qc_report}", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(bam)

    output:
    path "${id}_bam_stat.txt", emit: stats
    path "${id}_bam_flagstat.tsv", emit: flagstat

    script:
    """
    samtools stats ${bam} > ${id}_bam_stat.txt
    samtools flagstat ${bam} > ${id}_bam_flagstat.tsv
    """
}

//
process HAPLOTYPECALLER {
    publishDir "${params.outdir}/${params.vcf_hc}", mode: 'copy', overwrite: false

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
process DEEPVARIANT {
    publishDir "${params.outdir}/${params.vcf_dv}", mode: 'copy', overwrite: false

    scratch true

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id),
	path("${id}_dv.vcf.gz"),
	path("${id}_dv.vcf.gz.tbi"),
	emit: dv_vcf
    tuple val(id),
        path("${id}_dv.g.vcf.gz"),
        path("${id}_dv.g.vcf.gz.tbi"),
        emit: dv_gvcf

    script:
    """
    singularity exec \
      ${params.deepvariant} \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${params.reference_files.ref} \
        --reads=${bam} \
        --num_shards=${task.cpus} \
        --output_vcf=${id}_dv.vcf.gz \
        --output_gvcf=${id}_dv.g.vcf.gz
    """
}


//PANGENOME


process BUILD_REF_PATHS {

    input:
    path gbz

    output:
    path "ref_paths.txt", emit: ref_paths

    script:
    """
    # Extract paths and sort in one pipe to avoid creating large intermediate files
    vg paths -L -x ${gbz} | grep '^GRCh38' > ref_paths.txt
    """
}


//
process GIRAFFE {
    scratch true   // ⭐ensures node-local temp cleanup

    input:
    tuple val(id), val(lane), path(r1), path(r2)
    path ref_paths

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_p_markdup.bam"),
          path("${id}_${lane}_p_markdup.bam.bai"),
          emit: giraffe_bam

    script:
    """
    set -euo pipefail

    # 1. Giraffe alignment (BAM output)
    vg giraffe -t ${task.cpus} \
       -Z ${params.reference_pangenome.gbz} \
       -d ${params.reference_pangenome.dist} \
       -z ${params.reference_pangenome.zipcodes} \
       -m ${params.reference_pangenome.min} \
       -f ${r1} -f ${r2} \
       --ref-paths ${ref_paths} \
       --max-tail-length 50 \
       --max-chain-connection 50 \
       --wfa-max-mismatches 2 \
       --output-format BAM | \
    samtools sort -n -@ ${task.cpus} -m 4G -u -T ${params.tmp_dir}/${id}_${lane}_pnsort - | \
    samtools fixmate -m -u -@ ${task.cpus} - - | \
    samtools sort -@ ${task.cpus} -m 4G -T ${params.tmp_dir}/${id}_${lane}_pcsort - | \
    samtools markdup -@ ${task.cpus} - - | \
    samtools view -h - | grep -v "CHM13#0#" | \
    samtools view -b | \
    samtools reheader -c 'sed "s/GRCh38#0#//g"' - | \
    samtools addreplacerg \
       -r '@RG\\tID:${id}.${lane}\\tSM:${id}\\tLB:lib1\\tPL:ILLUMINA\\tPU:${lane}' \
       -o ${id}_${lane}_p_markdup.bam -

    # 2. Index final BAM
    samtools index ${id}_${lane}_p_markdup.bam
    """
}


//
process pBASERECALIBRATOR {
    scratch true   // ⭐ensures node-local temp cleanup
    
    input:
    tuple val(id), val(lane), path(bam), path(bai)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_p.table"), 
          emit: p_table

    script:
    """
    gatk BaseRecalibrator \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        --known-sites ${params.reference_files.dbsnp} \
        -O ${id}_${lane}_p.table
    """
}


//
process pApplyBQSR {
    input:
    tuple val(id), val(lane), path(bam), path(bai), path(table)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_p_bqsr.bam"),
          path("${id}_${lane}_p_bqsr.bam.bai"),
          emit: p_bqsr_bam

    script:
    """
    gatk ApplyBQSR \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        --bqsr-recal-file ${table} \
        -O ${id}_${lane}_p_bqsr.bam

    samtools index ${id}_${lane}_p_bqsr.bam
    """
}

//
process pMERGE_BAMS {
    publishDir "${params.outdir}/${params.bam}", mode: 'copy', overwrite: false
    scratch true

    input:
    tuple val(id), path(bams)

    output:
    tuple val(id),
          path("${id}_p_merged.bam"),
          path("${id}_p_merged.bam.bai"),
          emit: p_merged_bam

    script:
    """
    samtools merge -@ ${task.cpus} \
          -o ${id}_p_merged.bam \
             ${bams.join(' ')}
    samtools index ${id}_p_merged.bam
    """
}

//
process pBAM_QC {
    publishDir "${params.outdir}/${params.bam_qc_report}", mode: 'copy', overwrite: false

    input:
    tuple val(id), path(bam)

    output:
    path "${id}_p_bam_stat.txt", emit: stats
    path "${id}_p_bam_flagstat.tsv", emit: flagstat

    script:
    """
    samtools stats ${bam} > ${id}_p_bam_stat.txt
    samtools flagstat ${bam} > ${id}_p_bam_flagstat.tsv
    """
}



//
process pHAPLOTYPECALLER {
    publishDir "${params.outdir}/${params.vcf_p_hc}", mode: 'copy', overwrite: false
    scratch true

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id),
        path("${id}_hc_p.g.vcf.gz"),
        path("${id}_hc_p.g.vcf.gz.tbi"),
        emit: hc_p_gvcf

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus} \
        -O ${id}_hc_p.g.vcf.gz
    """
}

//
process pDEEPVARIANT {
    publishDir "${params.outdir}/${params.vcf_p_dv}", mode: 'copy', overwrite: false
    scratch true

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id),
        path("${id}_dv_p.vcf.gz"),
        path("${id}_dv_p.vcf.gz.tbi"),
        emit: dv_p_vcf
    tuple val(id),          
        path("${id}_dv_p.g.vcf.gz"),
        path("${id}_dv_p.g.vcf.gz.tbi"),
        emit: dv_p_gvcf

    script:
    """
    singularity exec \
      ${params.deepvariant} \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${params.reference_files.ref} \
        --reads=${bam} \
        --num_shards=${task.cpus} \
        --output_vcf=${id}_dv_p.vcf.gz \
        --output_gvcf=${id}_dv_p.g.vcf.gz
    """
}

//
workflow {
    // 1. Setup Input Channels
    ch_samples = Channel.fromPath(params.input).splitCsv(header: true).map { row -> tuple(row.sample, row.lane, file(row.fastq_1), file(row.fastq_2)) }
    
    //
    ch_chroms = Channel.fromList((1..22).collect{"chr$it"} + ['chrX','chrY','chrM'])
    
    // 2. Linear Track Pre-processing
    fastqc_ch    = FASTQC(ch_samples)
    //trimmed_ch = TRIMGALORE(ch_samples).trimmed_fq

    bwa_ch    = BWA(ch_samples).markdup_bam
    tables_ch = BASERECALIBRATOR(bwa_ch).table
    
    // Join original BAM with the specific chromosome table
    apply_input   = bwa_ch.join(tables_ch, by: [0,1])

    ch_bqsr_bams  = ApplyBQSR(apply_input).bqsr_bam
    
    // Merge scattered BAMs by Sample ID
    ch_merged_in  = ch_bqsr_bams.map{ id, lane, bam, bai -> [id, bam] }.groupTuple()
    merged_bam_ch = MERGE_BAMS(ch_merged_in).merged_bam
    ch_bam_qc     = BAM_QC(merged_bam_ch)


    // 3. Linear Variant Calling
    // Use .combine() instead of .cross() for more reliable scatter/parallelism
    hc_varcall = HAPLOTYPECALLER(merged_bam_ch)
    dv_varcall = DEEPVARIANT(merged_bam_ch)

    // 4. Pangenome Track
    ref_paths_ch = BUILD_REF_PATHS(file(params.reference_pangenome.gbz)).ref_paths

    giraffe_ch   = GIRAFFE(ch_samples, ref_paths_ch).giraffe_bam

    p_tables_ch = pBASERECALIBRATOR(giraffe_ch).p_table

    // Join original BAM with the specific id and lane table
    p_apply_input  = giraffe_ch.join(p_tables_ch, by: [0,1])
    
    ch_bqsr_bams_p = pApplyBQSR(p_apply_input).p_bqsr_bam

    // Merge scattered BAMs by Sample ID
    p_ch_merged_in  = ch_bqsr_bams_p.map{ id, lane, bam, bai -> [id, bam] }.groupTuple()
    p_merged_bam_ch = pMERGE_BAMS(p_ch_merged_in).p_merged_bam
    p_ch_bam_qc     = pBAM_QC(p_merged_bam_ch)

    // 5. Pangenome Variant Calling
    p_hc_varcall = pHAPLOTYPECALLER(p_merged_bam_ch)
    p_dv_varcall = pDEEPVARIANT(p_merged_bam_ch)
}
