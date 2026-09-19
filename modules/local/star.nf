process STAR_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/star", mode: params.publish_mode

    input:
    tuple val(meta), path(prepare_dir)
    path mapper_index

    output:
    tuple val(meta), path('nelrune.mapper.bam'), emit: bam
    tuple val(meta), path('Log.final.out'), path('Log.out'), path('Log.progress.out'), emit: qc

    script:
    def starBin = params.mapper_bin ?: 'STAR'
    def mapperOptions = params.mapper_options ?: ''

    """
    mapfile -t prepared_fastqs < <(find ${prepare_dir}/prepared_fastqs -maxdepth 1 -type f -name 'prepared.thread-*.fastq.gz' -print | sort)
    if [[ \${#prepared_fastqs[@]} -eq 0 ]]; then
        echo 'NELRUNE_PREPARE produced no prepared FASTQ shards' >&2
        exit 1
    fi
    read_files=\$(IFS=,; echo "\${prepared_fastqs[*]}")

    ${starBin} \\
        --runThreadN ${params.mapper_threads} \\
        --genomeDir ${mapper_index} \\
        --readFilesIn "\$read_files" \\
        --readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes NH HI AS nM \\
        --outReadsUnmapped None \\
        --outFileNamePrefix star. \\
        ${mapperOptions}

    mv star.Aligned.out.bam nelrune.mapper.bam
    mv star.Log.final.out Log.final.out
    mv star.Log.out Log.out
    mv star.Log.progress.out Log.progress.out
    """

    stub:
    """
    touch nelrune.mapper.bam Log.final.out Log.out Log.progress.out
    """
}
