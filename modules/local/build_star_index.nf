process BUILD_STAR_INDEX {
    tag "${genome.simpleName}"
    publishDir { index_publish_dir }, mode: params.reference_publish_mode

    input:
    path genome
    path gtf
    val index_publish_dir
    val index_name

    output:
    path "${index_name}", emit: index

    script:
    """
    mkdir -p ${index_name}

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index_name} \\
        --genomeFastaFiles ${genome} \\
        --sjdbGTFfile ${gtf}
    """
    stub:
    """
    mkdir -p ${index_name}
    touch ${index_name}/Genome
    touch ${index_name}/SA
    touch ${index_name}/SAindex
    touch ${index_name}/genomeParameters.txt
    """

}
