process BUILD_MINIMAP2_INDEX {
    tag "${genome.simpleName}"
    publishDir { index_publish_dir }, mode: params.reference_publish_mode
    input:
    path genome
    val index_publish_dir
    val index_name
    output:
    path "${index_name}", emit: index
    script:
    """
    minimap2 -t ${task.cpus} -x splice -d ${index_name} ${genome}
    """
    stub:
    """
    printf 'stub minimap2 index\n' > ${index_name}
    """
}
