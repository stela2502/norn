process NELRUNE {
    tag "${meta.id}"
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
    """
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
        --outpath nelrune_out \\
        --no-health-server
    """
}
