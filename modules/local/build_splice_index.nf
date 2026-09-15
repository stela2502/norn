process BUILD_SPLICE_INDEX {
    tag "${gtf.simpleName}"
    publishDir { index_publish_dir }, mode: params.reference_publish_mode

    input:
    path gtf
    val index_publish_dir
    val index_name
    val lumrik_runtime

    output:
    path "${index_name}", emit: index
    path "${index_name}.lumrik-runtime", emit: runtime_marker

    script:
    """
    gtf-splice-index build \\
        --annotation ${gtf} \\
        --index ${index_name}
    printf '%s\\n' '${lumrik_runtime}' > ${index_name}.lumrik-runtime
    """
    stub:
    """
    touch ${index_name}
    printf '%s\\n' '${lumrik_runtime}' > ${index_name}.lumrik-runtime
    """

}
