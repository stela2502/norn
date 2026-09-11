process CREATE_SEURAT {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/seurat", mode: params.publish_mode

    container params.seurat_container

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta), path("${meta.id}.seurat.rds"), emit: object

    script:
    """
    vdj_args=()
    if [[ -f vdj_calls.tsv && -f vdj_receptors.tsv ]]; then
        vdj_args+=(--vdj-calls vdj_calls.tsv --vdj-receptors vdj_receptors.tsv)
    fi

    Rscript \$(command -v create_seurat.R) \
        --exonic exonic \
        --intronic intronic \
        --project '${meta.id}' \
        --output '${meta.id}.seurat.rds' \
        "\${vdj_args[@]}"
    """
    stub:
    """
    touch '${meta.id}.seurat.rds'
    """

}
