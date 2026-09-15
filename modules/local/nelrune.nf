process NELRUNE {
    tag "${meta.id}"
    debug { params.health_server as boolean }
    publishDir "${params.outdir}/${meta.id}/nelrune", mode: params.publish_mode

    input:
    tuple val(meta), path(r1), path(r2)
    path splice_index
    path mapper_index

    output:
    tuple val(meta), path('nelrune_out/exonic'), emit: exonic
    tuple val(meta), path('nelrune_out/intronic'), emit: intronic
    tuple val(meta), path('nelrune_out/nelrune.mapper.bam'), emit: bam
    tuple val(meta), path('nelrune_out/nelrune-report.txt'), path('nelrune_out/nelrune.log'), path('nelrune_out/nelrune.metrics.tsv'), emit: qc

    script:
    def healthPort = (params.nelrune_health_port_base as int) + (task.index as int) - 1
    def healthArg = params.health_server ? "--health-port ${healthPort} --health-hostname \"\\$health_host\"" : '--no-health-server'
    def arg = { flag, value -> value != null && value.toString() != '' ? "${flag} ${value}" : '' }
    def boolArg = { flag, value -> value as boolean ? flag : '' }
    def listArg = { flag, value ->
        if (value == null || value.toString().trim() == '') return ''
        def vals = value instanceof Collection ? value : value.toString().split(',').collect { it.trim() }.findAll { it }
        vals ? "${flag} ${vals.join(' ')}" : ''
    }

    def primerStructureArg = arg('--primer-structure', params.primer_structure)
    def whitelistArg = arg('--whitelist', params.whitelist)
    def mapperBinArg = arg('--mapper-bin', params.mapper_bin)
    def mapperOptionsArg = params.mapper_options != null && params.mapper_options.toString() != '' ? "--mapper-options '${params.mapper_options.toString().replace("'", "'\\''")}'" : ''
    def genomeArg = params.genome ? "--genome ${params.genome}" : ''
    def vcfArg = arg('--vcf', params.vcf)
    def maxReadsArg = arg('--max-reads', params.max_reads)
    def readTagArg = listArg('--read-tag-table', params.read_tag_table)
    def bamOutArg = arg('--bam-out', params.bam_out)
    def featuresArg = listArg('--additional-features', params.additional_features)

    """
    mkdir -p lumrik_tmp
    export LUMRIK_TMPDIR="\$PWD/lumrik_tmp"
    health_host="\$(hostname -s 2>/dev/null || hostname)"

    nelrune \\
        --r1 ${r1.join(' ')} \\
        --r2 ${r2.join(' ')} \\
        --chemistry ${meta.chemistry} \\
        ${primerStructureArg} \\
        ${whitelistArg} \\
        --whitelist-mismatches ${params.whitelist_mismatches} \\
        ${boolArg('--detect-reverse-complement', params.detect_reverse_complement)} \\
        --mapper ${params.mapper} \\
        ${mapperBinArg} \\
        --mapper-index ${mapper_index} \\
        ${mapperOptionsArg} \\
        --mapper-threads ${params.mapper_threads} \\
        ${boolArg('--mapper-paired', params.mapper_paired)} \\
        ${boolArg('--mapper-keep-multimappers', params.mapper_keep_multimappers)} \\
        --index ${splice_index} \\
        ${genomeArg} \\
        ${vcfArg} \\
        --quant-mode ${params.quant_mode} \\
        --min-mapq ${params.min_mapq} \\
        ${maxReadsArg} \\
        ${boolArg('--read1-only', params.read1_only)} \\
        ${boolArg('--no-genome-refine', params.no_genome_refine)} \\
        ${boolArg('--require-strand', params.require_strand)} \\
        ${boolArg('--require-exact-junction-chain', params.require_exact_junction_chain)} \\
        --max-5p-overhang-bp ${params.max_5p_overhang_bp} \\
        --max-3p-overhang-bp ${params.max_3p_overhang_bp} \\
        --allowed-intronic-gap-size ${params.allowed_intronic_gap_size} \\
        --snp-min-anchor ${params.snp_min_anchor} \\
        ${readTagArg} \\
        --rt-read-id-column ${params.rt_read_id_column} \\
        --rt-cell-column ${params.rt_cell_column} \\
        --rt-cell-qual-column ${params.rt_cell_qual_column} \\
        --rt-umi-column ${params.rt_umi_column} \\
        --rt-umi-qual-column ${params.rt_umi_qual_column} \\
        --rt-original-read-id-column ${params.rt_original_read_id_column} \\
        --analysis-type ${params.analysis_type} \\
        --cell-tag ${params.cell_tag} \\
        --umi-tag ${params.umi_tag} \\
        ${bamOutArg} \\
        ${featuresArg} \\
        --additional-feature-min-hits ${params.additional_feature_min_hits} \\
        --min-insert-len ${params.min_insert_len} \\
        --min-transcript-len ${params.min_transcript_len} \\
        --min-cell-counts ${params.min_cell_counts} \\
        --threads ${params.nelrune_threads} \\
        --outpath nelrune_out \\
        ${healthArg} \\
        ${params.nelrune_extra_args ?: ''} \\
        >nelrune.console.log 2>&1 &
    lumrik_pid=\$!

    if ${params.health_server}; then
        health_ready=0
        for _ in \$(seq 1 150); do
            if ! kill -0 "\$lumrik_pid" 2>/dev/null; then break; fi
            if (echo > /dev/tcp/127.0.0.1/${healthPort}) >/dev/null 2>&1; then
                printf '\nNorn health server\n  NELRUNE (%s)  http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
                health_ready=1
                break
            fi
            sleep 0.2
        done
        if [[ "\$health_ready" -eq 0 ]] && kill -0 "\$lumrik_pid" 2>/dev/null; then
            printf '\nNorn health server\n  NELRUNE (%s)  starting on http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
        fi
    fi

    set +e
    wait "\$lumrik_pid"
    lumrik_status=\$?
    set -e
    if [[ "\$lumrik_status" -ne 0 ]]; then
        tail -n 80 nelrune.console.log >&2 || true
        exit "\$lumrik_status"
    fi
    """
    stub:
    """
    mkdir -p nelrune_out/exonic nelrune_out/intronic
    touch nelrune_out/exonic/matrix.mtx.gz nelrune_out/exonic/features.tsv.gz nelrune_out/exonic/barcodes.tsv.gz
    touch nelrune_out/intronic/matrix.mtx.gz nelrune_out/intronic/features.tsv.gz nelrune_out/intronic/barcodes.tsv.gz
    touch nelrune_out/nelrune.mapper.bam
    touch nelrune_out/nelrune-report.txt nelrune_out/nelrune.log nelrune_out/nelrune.metrics.tsv
    """
}
