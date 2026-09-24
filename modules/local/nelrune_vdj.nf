process NELRUNE_VDJ {
    tag "${meta.id}"
    debug { params.health_server as boolean }
    publishDir "${params.outdir}/${meta.id}/vdj", mode: params.publish_mode

    input:
    tuple val(meta), path(bam)
    path vdj_index

    output:
    tuple val(meta), path('vdj_out/vdj_calls.tsv'), path('vdj_out/vdj_receptors.tsv'), path('vdj_out/airr_rearrangements.tsv'), path('vdj_out/vdj-mapping-info.txt'), emit: tables
    tuple val(meta), path('vdj_out/vdj_observed.fasta'), path('vdj_out/vdj_naive.fasta'), optional: true, emit: sequences

    script:
    def bdVersions = [
        'bd-v1': 'v1',
        'bd-v2-96': 'v2.96',
        'bd-v2-384': 'v2.384'
    ]
    def bdVersion = params.vdj_bd_cell_version ?: bdVersions[meta.chemistry]
    def bdArg = bdVersion ? "--bd-cell-version ${bdVersion}" : ''
    def seqArg = params.vdj_write_sequences ? '--write-sequences' : ''
    def gtfArg = params.gtf ? "--gtf ${params.gtf}" : ''
    def genomeArg = params.genome ? "--genome ${params.genome}" : ''
    def healthPort = (params.vdj_health_port_base as int) + (task.index as int) - 1
    def healthArg = params.health_server ? "--health-port ${healthPort} --health-hostname \"\$health_host\"" : '--no-health-server'
    """
    health_host="\$(hostname -s 2>/dev/null || hostname)"

    nelrune-vdj \\
        --bam ${bam} \\
        --index ${vdj_index} \\
        --out vdj_out \\
        --threads ${params.vdj_threads} \\
        ${bdArg} \\
        ${seqArg} \\
        ${healthArg} \\
        >nelrune-vdj.console.log 2>&1 &
    lumrik_pid=\$!

    if ${params.health_server}; then
        health_ready=0
        for _ in \$(seq 1 150); do
            if ! kill -0 "\$lumrik_pid" 2>/dev/null; then
                break
            fi
            if (echo > /dev/tcp/127.0.0.1/${healthPort}) >/dev/null 2>&1; then
                printf '\nNorn health server\n  NELRUNE_VDJ (%s)  http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
                health_ready=1
                break
            fi
            sleep 0.2
        done
        if [[ "\$health_ready" -eq 0 ]] && kill -0 "\$lumrik_pid" 2>/dev/null; then
            printf '\nNorn health server\n  NELRUNE_VDJ (%s)  starting on http://%s:%s\n\n' '${meta.id}' "\$health_host" '${healthPort}'
        fi
    fi

    set +e
    wait "\$lumrik_pid"
    lumrik_status=\$?
    set -e
    if [[ "\$lumrik_status" -ne 0 ]]; then
        tail -n 80 nelrune-vdj.console.log >&2 || true
        exit "\$lumrik_status"
    fi
    """
    stub:
    """
    mkdir -p vdj_out
    printf 'cell\trecombination_id\tproductivity_status\n' > vdj_out/vdj_calls.tsv
    printf 'cell\theavy_recombination_id\tlight_recombination_id\n' > vdj_out/vdj_receptors.tsv
    touch vdj_out/airr_rearrangements.tsv vdj_out/vdj-mapping-info.txt
    """

}
