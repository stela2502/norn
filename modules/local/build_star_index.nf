process BUILD_STAR_INDEX {
    tag "${genome.simpleName}"
    publishDir "${params.outdir}/reference/star", mode: params.reference_publish_mode

    input:
    path genome
    path gtf

    output:
    path 'star_index', emit: index

    script:
    """
    mkdir -p star_index

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${genome} \\
        --sjdbGTFfile ${gtf}
    """
}
