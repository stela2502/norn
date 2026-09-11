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
    """
    nelrune-vdj \\
        --exonic ${exonic} \\
        --bam ${bam} \\
        --index ${vdj_index} \\
        --out vdj_out \\
        --threads ${params.vdj_threads} \\
        ${bdArg} \\
        ${seqArg} \\
        --no-health-server
    """
}
