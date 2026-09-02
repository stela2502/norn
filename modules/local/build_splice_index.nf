process BUILD_SPLICE_INDEX {
    tag "${gtf.simpleName}"
    publishDir "${params.outdir}/reference/splice", mode: params.reference_publish_mode

    input:
    path gtf

    output:
    path 'reference.splice.idx', emit: index

    script:
    """
    gtf-splice-index build \\
        --annotation ${gtf} \\
        --index reference.splice.idx
    """
}
