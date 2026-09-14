process BUILD_VDJ_INDEX {
    tag "${gtf.simpleName}"
    publishDir { index_publish_dir }, mode: params.reference_publish_mode

    input:
    path gtf
    path genome
    val index_publish_dir
    val index_name

    output:
    path "${index_name}", emit: index

    script:
    """
    vdj-index \\
        --gtf ${gtf} \\
        --genome ${genome} \\
        --out ${index_name}
    """
    stub:
    """
    touch ${index_name}
    """

}
