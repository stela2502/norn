process BUILD_SPLICE_INDEX {
    tag "${gtf.simpleName}"
    publishDir { index_publish_dir }, mode: params.reference_publish_mode

    input:
    path gtf
    val index_publish_dir
    val index_name

    output:
    path "${index_name}", emit: index

    script:
    """
    gtf-splice-index build \\
        --annotation ${gtf} \\
        --index ${index_name}
    """
    stub:
    """
    touch ${index_name}
    """

}
