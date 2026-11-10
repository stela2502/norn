process NELRUNE_ONT {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/nelrune", mode: params.publish_mode
    input:
    tuple val(meta), path(bam)
    path mapper_index
    path splice_index
    path genome
    path vcf
    output:
    tuple val(meta), path('nelrune_out/filtered/exprs'), emit: filtered_exprs
    // Publish the complete Nelrune output, including novel transcripts and SNP evidence.
    tuple val(meta), path('nelrune_out'), emit: full_results
    tuple val(meta), path('nelrune_out/nelrune-report.txt'), optional: true, emit: qc
    script:
    def opt = { flag, value -> value != null && value.toString().trim() ? "${flag} ${value}" : '' }
    def flag = { name, value -> value as boolean ? name : '' }
    def vcfArg = vcf ? "--vcf ${vcf}" : ''
    def healthPort = (params.nelrune_health_port_base as int) + (task.index as int) - 1
    def healthArg = params.health_server ? "--health-port ${healthPort}" : '--no-health-server'
    """
    mkdir -p lumrik_tmp
    export LUMRIK_TMPDIR="\$PWD/lumrik_tmp"
    nelrune \
        --bam ${bam} \
        --chemistry ${meta.chemistry} \
        --mapper ${params.mapper} \
        --mapper-index ${mapper_index} \
        --mapper-threads ${params.mapper_threads} \
        --threads ${params.nelrune_quant_threads} \
        --index ${splice_index} \
        --genome ${genome} \
        --mapper-keep-multimappers \
        ${opt('--mapper-options', params.mapper_options)} \
        ${opt('--primer-structure', params.primer_structure)} \
        ${opt('--whitelist', params.whitelist)} \
        --whitelist-mismatches ${params.whitelist_mismatches} \
        ${flag('--detect-reverse-complement', params.detect_reverse_complement)} \
        ${flag('--barcode-crosscheck', params.barcode_crosscheck)} \
        ${flag('--ont-diagnostics', params.ont_diagnostics)} \
        ${opt('--dump-failed-reads', params.dump_failed_reads)} \
        ${vcfArg} \
        ${opt('--max-reads', params.max_reads)} \
        ${opt('--min-umi-count', params.min_umi_counts)} \
        --quant-mode ${params.quant_mode} \
        --min-mapq ${params.min_mapq} \
        --snp-min-anchor ${params.snp_min_anchor} \
        ${healthArg} \
        ${params.nelrune_ont_extra_args ?: ''} \
        --outpath nelrune_out
    """
    stub:
    """
    mkdir -p nelrune_out/filtered/exprs
    touch nelrune_out/filtered/exprs/matrix.mtx.gz nelrune_out/filtered/exprs/barcodes.tsv.gz nelrune_out/filtered/exprs/features.tsv.gz nelrune_out/nelrune-report.txt
    """
}
