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
    tuple val(id), val(lane), val(chrom), path(r1), path(r2)

    output:
    tuple val(id), val(lane), val(chrom),
          path("${id}_${lane}_${chrom}_markdup.bam"),
          path("${id}_${lane}_${chrom}_markdup.bam.bai"),
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
       -o ${id}_${lane}_${chrom}_markdup.bam -
	
    # 2. Indexing BAM file
    samtools index ${id}_${lane}_${chrom}_markdup.bam
    """
}

//
process BASERECALIBRATOR {
    scratch true   // ⭐ensures node-local temp cleanup

    input:
    tuple val(id), val(lane), val(chrom), path(bam), path(bai)

    output:
    tuple val(id), val(lane), val(chrom),
          path("${id}_${lane}_${chrom}.table"), 
          emit: table
    
    script:
    """
    gatk BaseRecalibrator \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        --known-sites ${params.reference_files.dbsnp} \
        -L ${chrom} \
        -O ${id}_${lane}_${chrom}.table
    """
}

//
process ApplyBQSR {
    input:
    tuple val(id), val(lane), val(chrom), path(bam), path(bai), path(table)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_${chrom}_bqsr.bam"),
          path("${id}_${lane}_${chrom}_bqsr.bam.bai"),
          emit: bqsr_bam

    script:
    """
    gatk ApplyBQSR \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -L ${chrom} \
        --bqsr-recal-file ${table} \
        -O ${id}_${lane}_${chrom}_bqsr.bam

    samtools index ${id}_${lane}_${chrom}_bqsr.bam
    """
}


//
process MERGE_BAMS {
    publishDir "${params.outdir}/bams", mode: 'copy', overwrite: false
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
process HAPLOTYPECALLER {
    scratch true

    input:
    tuple val(id), path(bam), path(bai), val(chrom)

    output:
    tuple val(id),
        path("${id}_${chrom}_hc.g.vcf.gz"),
        path("${id}_${chrom}_hc.g.vcf.gz.tbi"),
        emit: hc_gvcf
          
    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -L ${chrom} \
        -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus} \
        -O ${id}_${chrom}_hc.g.vcf.gz
    """
}

//
process DEEPVARIANT {
    scratch true

    input:
    tuple val(id), path(bam), path(bai), val(chrom)

    output:
    tuple val(id),
	path("${id}_${chrom}_dv.vcf.gz"),
	path("${id}_${chrom}_dv.vcf.gz.tbi"),
	emit: dv_vcf
    tuple val(id),
        path("${id}_${chrom}_dv.g.vcf.gz"),
        path("${id}_${chrom}_dv.g.vcf.gz.tbi"),
        emit: dv_gvcf

    script:
    """
    singularity exec \
      ${params.deepvariant} \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${params.reference_files.ref} \
        --reads=${bam} \
        --regions=${chrom} \
        --num_shards=${task.cpus} \
        --output_vcf=${id}_${chrom}_dv.vcf.gz \
        --output_gvcf=${id}_${chrom}_dv.g.vcf.gz
    """
}

//
process CONCAT_HC_GVCFS {
    publishDir "${params.outdir}/hc_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(gvcfs), path(tbis)

    output:
    tuple val(id),
       path("${id}_hc_concat.g.vcf.gz"),
       path("${id}_hc_concat.g.vcf.gz.tbi"),
       emit: hc_gvcf
       path("hc_files.list")

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${gvcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> hc_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f hc_files.list -O z -o ${id}_hc_concat.g.vcf.gz
    bcftools index -t ${id}_hc_concat.g.vcf.gz
    """
}

//
process CONCAT_DV_VCFS {
    publishDir "${params.outdir}/dv_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(vcfs), path(tbis)

    output:
    tuple val(id),
        path("${id}_dv_concat.vcf.gz"),
        path("${id}_dv_concat.vcf.gz.tbi"),
        emit: dv_vcf
        path("dv_files.list")

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${vcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> dv_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f dv_files.list -O z -o ${id}_dv_concat.vcf.gz
    bcftools index -t ${id}_dv_concat.vcf.gz
    """
}


//
process CONCAT_DV_GVCFS {
    publishDir "${params.outdir}/dv_g_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(gvcfs), path(tbis)

    output:
    tuple val(id),
        path("${id}_dv_concat.g.vcf.gz"),
        path("${id}_dv_concat.g.vcf.gz.tbi"),
        emit: dv_gvcf
        path("dv_g_files.list")

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${gvcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> dv_g_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f dv_g_files.list -O z -o ${id}_dv_concat.g.vcf.gz
    bcftools index -t ${id}_dv_concat.g.vcf.gz
    """
}



// PANGENOME



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
    tuple val(id), val(lane), val(chrom), path(r1), path(r2)
    path ref_paths

    output:
    tuple val(id), val(lane), val(chrom),
          path("${id}_${lane}_${chrom}_p_markdup.bam"),
          path("${id}_${lane}_${chrom}_p_markdup.bam.bai"),
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
       -o ${id}_${lane}_${chrom}_p_markdup.bam -

    # 2. Index final BAM
    samtools index ${id}_${lane}_${chrom}_p_markdup.bam
    """
}


//
process pBASERECALIBRATOR {
    scratch true   // ⭐ensures node-local temp cleanup
    
    input:
    tuple val(id), val(lane), val(chrom), path(bam), path(bai)

    output:
    tuple val(id), val(lane), val(chrom),
          path("${id}_${lane}_${chrom}_p.table"), 
          emit: p_table

    script:
    """
    gatk BaseRecalibrator \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -L ${chrom} \
        --known-sites ${params.reference_files.dbsnp} \
        -O ${id}_${lane}_${chrom}_p.table
    """
}


//
process pApplyBQSR {
    input:
    tuple val(id), val(lane), val(chrom), path(bam), path(bai), path(table)

    output:
    tuple val(id), val(lane),
          path("${id}_${lane}_${chrom}_p_bqsr.bam"),
          path("${id}_${lane}_${chrom}_p_bqsr.bam.bai"),
          emit: p_bqsr_bam

    script:
    """
    gatk ApplyBQSR \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -L ${chrom} \
        --bqsr-recal-file ${table} \
        -O ${id}_${lane}_${chrom}_p_bqsr.bam

    samtools index ${id}_${lane}_${chrom}_p_bqsr.bam
    """
}

//
process pMERGE_BAMS {
    publishDir "${params.outdir}/bams", mode: 'copy', overwrite: false
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
process pHAPLOTYPECALLER {
    scratch true

    input:
    tuple val(id), path(bam), path(bai), val(chrom)

    output:
    tuple val(id),
        path("${id}_${chrom}_hc_p.g.vcf.gz"),
        path("${id}_${chrom}_hc_p.g.vcf.gz.tbi"),
        emit: hc_p_gvcf

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference_files.ref} \
        -I ${bam} \
        -L ${chrom} \
        -ERC GVCF \
        --native-pair-hmm-threads ${task.cpus} \
        -O ${id}_${chrom}_hc_p.g.vcf.gz
    """
}

//
process pDEEPVARIANT {
    scratch true

    input:
    tuple val(id), path(bam), path(bai), val(chrom)

    output:
    tuple val(id),
        path("${id}_${chrom}_dv_p.vcf.gz"),
        path("${id}_${chrom}_dv_p.vcf.gz.tbi"),
        emit: dv_p_vcf
    tuple val(id),          
        path("${id}_${chrom}_dv_p.g.vcf.gz"),
        path("${id}_${chrom}_dv_p.g.vcf.gz.tbi"),
        emit: dv_p_gvcf

    script:
    """
    singularity exec \
      ${params.deepvariant} \
      /opt/deepvariant/bin/run_deepvariant \
        --model_type=WGS \
        --ref=${params.reference_files.ref} \
        --reads=${bam} \
        --regions=${chrom} \
        --num_shards=${task.cpus} \
        --output_vcf=${id}_${chrom}_dv_p.vcf.gz \
        --output_gvcf=${id}_${chrom}_dv_p.g.vcf.gz
    """
}

//
process pCONCAT_HC_GVCFS {
    publishDir "${params.outdir}/hc_p_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(gvcfs), path(tbis)

    output:
    tuple val(id),
       path("${id}_hc_p_concat.g.vcf.gz"),
       path("${id}_hc_p_concat.g.vcf.gz.tbi"),
       emit: hc_p_gvcf

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${gvcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> hc_p_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f hc_p_files.list -O z -o ${id}_hc_p_concat.g.vcf.gz
    bcftools index -t ${id}_hc_p_concat.g.vcf.gz
    """
}

//
process pCONCAT_DV_VCFS {
    publishDir "${params.outdir}/dv_p_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(vcfs), path(tbis)

    output:
    tuple val(id),
        path("${id}_dv_p_concat.vcf.gz"),
        path("${id}_dv_p_concat.vcf.gz.tbi"),
        emit: dv_p_vcf
        path("dv_p_files.list")

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${vcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> dv_p_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f dv_p_files.list -O z -o ${id}_dv_p_concat.vcf.gz
    bcftools index -t ${id}_dv_p_concat.vcf.gz
    """
}

//
process pCONCAT_DV_GVCFS {
    //publishDir "${params.outdir}/dv_p_g_concat", mode: 'copy', overwrite: false
    input:
    tuple val(id), path(gvcfs), path(tbis)

    output:
    tuple val(id),
        path("${id}_dv_p_concat.g.vcf.gz"),
        path("${id}_dv_p_concat.g.vcf.gz.tbi"),
        emit: dv_p_gvcf
        path("dv_p_g_files.list")

    script:
    """
    set -euo pipefail
    # 1. Define the desired order
    # This creates a list: chr1, chr2 ... chr22, chrX, chrY, chrM
    order=\$(printf "chr%s\\n" {1..22} X Y M)

    # 2. Match actual files against that order to create the file list
    for c in \$order; do
        for f in ${gvcfs.join(' ')}; do
            if [[ "\$f" == *"\${c}_"* ]]; then
                echo "\$f" >> dv_p_g_files.list
                break
            fi
        done
    done

    # 3. Concatenate and index
    bcftools concat -f dv_p_g_files.list -O z -o ${id}_dv_p_concat.g.vcf.gz
    bcftools index -t ${id}_dv_p_concat.g.vcf.gz
    """
}

//
workflow {
    // 1. Setup Input Channels
    ch_samples = Channel.fromPath(params.input).splitCsv(header: true).map { row -> tuple(row.sample, row.lane, file(row.fastq_1), file(row.fastq_2)) }
    
    //
    ch_chroms = Channel.fromList((1..22).collect{"chr$it"} + ['chrX','chrY','chrM'])
    
    // 2. Linear Track Pre-processing
    //trimmed_ch = TRIMGALORE(ch_samples).trimmed_fq
    //bwa_ch     = BWA(trimmed_ch).markdup_bam
    sample_by_chrom = ch_samples.combine(ch_chroms).map { id, lane, r1, r2, chrom -> tuple(id, lane, chrom, r1, r2) }    
    bwa_ch = BWA(sample_by_chrom).markdup_bam

    // Combine aligned BAMs with chromosomes list to get new BAM
    //bwa_by_chrom = bwa_ch.combine(ch_chroms).map { id, lane, bam, bai, chrom -> tuple(id, lane, chrom, bam, bai) }

    //tables_ch   = BASERECALIBRATOR(bwa_by_chrom).table
    tables_ch   = BASERECALIBRATOR(bwa_ch).table
    
    // Join original BAM with the specific chromosome table
    //apply_input   = bwa_by_chrom.join(tables_ch, by: [0,1,2])
    apply_input   = bwa_ch.join(tables_ch, by: [0,1,2])

    ch_bqsr_bams  = ApplyBQSR(apply_input).bqsr_bam
    
    // Merge scattered BAMs by Sample ID
    ch_merged_in  = ch_bqsr_bams.map{ id, lane, bam, bai -> [id, bam] }.unique().groupTuple()
    merged_bam_ch = MERGE_BAMS(ch_merged_in).merged_bam

    // 3. Linear Variant Calling
    // Use .combine() instead of .cross() for more reliable scatter/parallelism
    hc_varcall = HAPLOTYPECALLER(merged_bam_ch.combine(ch_chroms))
    dv_varcall = DEEPVARIANT(merged_bam_ch.combine(ch_chroms))
    
    hc_gvcfs_ch = hc_varcall.hc_gvcf.map { id, gvcf, tbi -> tuple(id, gvcf, tbi) }.groupTuple()
    hc_concat   = CONCAT_HC_GVCFS(hc_gvcfs_ch)

    dv_vcfs_ch  = dv_varcall.dv_vcf.map { id, vcf, tbi -> tuple(id, vcf, tbi) }.groupTuple()
    dv_vcfs_concat = CONCAT_DV_VCFS(dv_vcfs_ch)

    dv_gvcfs_ch = dv_varcall.dv_gvcf.map { id, gvcf, tbi -> tuple(id, gvcf, tbi) }.groupTuple()
    dv_gvcf_concat = CONCAT_DV_GVCFS(dv_gvcfs_ch)


    // 4. Pangenome Track
    ref_paths_ch = BUILD_REF_PATHS(file(params.reference_pangenome.gbz)).ref_paths
    
    //giraffe_ch   = GIRAFFE(trimmed_ch, ref_paths_ch).giraffe_bam
    giraffe_ch   = GIRAFFE(sample_by_chrom, ref_paths_ch).giraffe_bam

    // Combine aligned BAMs with chromosomes list
    //giraffe_by_chrom = giraffe_ch.combine(ch_chroms).map { id, lane, bam, bai, chrom -> tuple(id, lane, chrom, bam, bai) }
    
    //p_tables_ch = pBASERECALIBRATOR(giraffe_by_chrom).p_table
    p_tables_ch = pBASERECALIBRATOR(giraffe_ch).p_table

    // Join original BAM with the specific chromosome table
    p_apply_input  = giraffe_ch.join(p_tables_ch, by: [0,1,2])
    
    ch_bqsr_bams_p = pApplyBQSR(p_apply_input).p_bqsr_bam

    // Merge scattered BAMs by Sample ID
    p_ch_merged_in  = ch_bqsr_bams_p.map{ id, lane, bam, bai -> [id, bam] }.unique().groupTuple()
    p_merged_bam_ch = pMERGE_BAMS(p_ch_merged_in).p_merged_bam

    // 5. Pangenome Variant Calling
    p_hc_varcall = pHAPLOTYPECALLER(p_merged_bam_ch.combine(ch_chroms))
    p_dv_varcall = pDEEPVARIANT(p_merged_bam_ch.combine(ch_chroms))

    p_hc_gvcfs_ch = p_hc_varcall.hc_p_gvcf.map { id, gvcf, tbi -> tuple(id, gvcf, tbi) }.groupTuple()
    p_hc_concat   = pCONCAT_HC_GVCFS(p_hc_gvcfs_ch)

    p_dv_vcfs_ch      = p_dv_varcall.dv_p_vcf.map { id, vcf, tbi -> tuple(id, vcf, tbi) }.groupTuple()
    p_dv_vcfs_concat  = pCONCAT_DV_VCFS(p_dv_vcfs_ch)

    p_dv_gvcfs_ch     = p_dv_varcall.dv_p_gvcf.map { id, gvcf, tbi -> tuple(id, gvcf, tbi) }.groupTuple()
    p_dv_gvcfs_concat = pCONCAT_DV_GVCFS(p_dv_gvcfs_ch)
}
