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
    def healthArg = params.health_server ? "--health-port ${healthPort} --health-hostname \"\$health_host\"" : '--no-health-server'
    def maxReadsArg = params.max_reads != null ? "--max-reads ${params.max_reads}" : ''
    """
    mkdir -p lumrik_tmp
    export LUMRIK_TMPDIR="\$PWD/lumrik_tmp"
    health_host="\$(hostname -s 2>/dev/null || hostname)"

    nelrune \\
        --r1 ${r1.join(' ')} \\
        --r2 ${r2.join(' ')} \\
        --chemistry ${meta.chemistry} \\
        --mapper ${params.mapper} \\
        --mapper-index ${mapper_index} \\
        --mapper-threads ${params.mapper_threads} \\
        --index ${splice_index} \\
        --threads ${params.nelrune_threads} \\
        --min-cell-counts ${params.min_cell_counts} \\
        ${maxReadsArg} \\
        --outpath nelrune_out \\
        ${healthArg} \\
        >nelrune.console.log 2>&1 &
    lumrik_pid=\$!

    if ${params.health_server}; then
        health_ready=0
        for _ in \$(seq 1 150); do
            if ! kill -0 "\$lumrik_pid" 2>/dev/null; then
                break
            fi
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
