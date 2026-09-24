process NELRUNE_PREPARE {
    tag "${meta.id}"
    debug { params.health_server as boolean }
    publishDir "${params.outdir}/${meta.id}/prepare", mode: params.publish_mode

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path('prepare_out'), emit: prepared

    script:
    def arg = { flag, value -> value != null && value.toString() != '' ? "${flag} ${value}" : '' }
    def boolArg = { flag, value -> value as boolean ? flag : '' }
    def listArg = { flag, value ->
        if (value == null || value.toString().trim() == '') return ''
        def vals = value instanceof Collection ? value : value.toString().split(',').collect { it.trim() }.findAll { it }
        vals ? "${flag} ${vals.join(' ')}" : ''
    }

    def primerStructureArg = arg('--primer-structure', params.primer_structure)
    def whitelistArg = arg('--whitelist', params.whitelist)
    def maxReadsArg = arg('--max-reads', params.max_reads)
    def featuresArg = listArg('--additional-features', params.additional_features)
    def healthPort = (params.nelrune_health_port_base as int) + (task.index as int) - 1
    def healthArg = params.health_server ? "--health-port ${healthPort} --health-hostname \"\$health_host\"" : '--no-health-server'

    """
    mkdir -p lumrik_tmp
    export LUMRIK_TMPDIR="\$PWD/lumrik_tmp"
    health_host="\$(hostname -s 2>/dev/null || hostname)"

    nelrune prepare-fastqs \\
        --r1 ${r1.join(' ')} \\
        --r2 ${r2.join(' ')} \\
        --chemistry ${meta.chemistry} \\
        ${primerStructureArg} \\
        ${whitelistArg} \\
        --whitelist-mismatches ${params.whitelist_mismatches} \\
        ${boolArg('--detect-reverse-complement', params.detect_reverse_complement)} \\
        ${featuresArg} \\
        --additional-feature-min-hits ${params.additional_feature_min_hits} \\
        --min-insert-len ${params.min_insert_len} \\
        ${maxReadsArg} \\
        --threads ${params.nelrune_prepare_threads} \\
        --outpath prepare_out \\
        ${healthArg} \\
        ${params.nelrune_prepare_extra_args ?: ''} \\
        >nelrune-prepare.console.log 2>&1 &
    lumrik_pid=\$!

    if ${params.health_server}; then
        health_ready=0
        for _ in \$(seq 1 150); do
            if ! kill -0 "\$lumrik_pid" 2>/dev/null; then break; fi
            if (echo > /dev/tcp/127.0.0.1/${healthPort}) >/dev/null 2>&1; then
                printf '\nNorn health server\n  NELRUNE_PREPARE (%s)  http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
                health_ready=1
                break
            fi
            sleep 0.2
        done
        if [[ "\$health_ready" -eq 0 ]] && kill -0 "\$lumrik_pid" 2>/dev/null; then
            printf '\nNorn health server\n  NELRUNE_PREPARE (%s)  starting on http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
        fi
    fi

    set +e
    wait "\$lumrik_pid"
    lumrik_status=\$?
    set -e
    if [[ "\$lumrik_status" -ne 0 ]]; then
        tail -n 80 nelrune-prepare.console.log >&2 || true
        exit "\$lumrik_status"
    fi
    """

    stub:
    """
    mkdir -p prepare_out/prepared_fastqs
    printf '@stub\\nACGT\\n+\\nIIII\\n' | gzip -c > prepare_out/prepared_fastqs/prepared.thread-000.fastq.gz
    printf 'format\\tnelrune-prepare-v1\\ncell_barcode_len\\t27\\nfastq\\tprepared_fastqs/prepared.thread-000.fastq.gz\\n' > prepare_out/prepare-manifest.tsv
    touch prepare_out/feature_observations.bin prepare_out/prepare-report.txt
    """
}
