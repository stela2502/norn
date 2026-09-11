process NELRUNE_VDJ {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/vdj", mode: params.publish_mode

    input:
    tuple val(meta), path(exonic), path(bam)
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
    def bdVersion = bdVersions[meta.chemistry]
    def bdArg = bdVersion ? "--bd-cell-version ${bdVersion}" : ''
    def seqArg = params.vdj_write_sequences ? '--write-sequences' : ''
    def healthPort = (params.vdj_health_port_base as int) + (task.index as int) - 1
    def healthArg = params.health_server ? "--health-port ${healthPort}" : '--no-health-server'
    """
    nelrune-vdj \\
        --exonic ${exonic} \\
        --bam ${bam} \\
        --index ${vdj_index} \\
        --out vdj_out \\
        --threads ${params.vdj_threads} \\
        ${bdArg} \\
        ${seqArg} \\
        ${healthArg}
    """
    stub:
    """
    mkdir -p vdj_out
    printf 'cell\trecombination_id\tproductivity_status\n' > vdj_out/vdj_calls.tsv
    printf 'cell\theavy_recombination_id\tlight_recombination_id\n' > vdj_out/vdj_receptors.tsv
    touch vdj_out/airr_rearrangements.tsv vdj_out/vdj-mapping-info.txt
    """

}
