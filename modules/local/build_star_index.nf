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

    star_genome="${genome}"
    case "${genome}" in
        *.gz)
            star_genome="star.genome.fa"
            gzip -dc "${genome}" > "\$star_genome"
            ;;
    esac

    star_gtf="${gtf}"
    case "${gtf}" in
        *.gz)
            star_gtf="star.annotation.gtf"
            gzip -dc "${gtf}" > "\$star_gtf"
            ;;
    esac

    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index_name} \\
        --genomeFastaFiles "\$star_genome" \\
        --sjdbGTFfile "\$star_gtf"
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
