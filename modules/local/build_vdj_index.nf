process BUILD_VDJ_INDEX {
    tag "${gtf.simpleName}"
    publishDir "${params.outdir}/reference/vdj", mode: params.reference_publish_mode

    input:
    path gtf
    path genome

    output:
    path 'reference.vdjidx', emit: index

    script:
    """
    vdj-index \\
        --gtf ${gtf} \\
        --genome ${genome} \\
        --out reference.vdjidx
    """
}
