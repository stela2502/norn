process CLONOMAP {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/clonomap", mode: params.publish_mode

    input:
    tuple val(meta), path(vdj_calls), path(vdj_receptors), path(airr), path(mapping_info)

    output:
    tuple val(meta), path('clonomap_out/cells.tsv'), emit: cells
    tuple val(meta), path('clonomap_out'), emit: results

    script:
    """
    def plots_arg = params.clonomap_plots ? '--plots' : ''

    clonomap_family --vdj-out . --out clonomap_out ${plots_arg} ${params.clonomap_extra_args ?: ''}
    """

    stub:
    """
    mkdir -p clonomap_out
    printf 'source\tcell\tfamily\tlc_clone\n' > clonomap_out/cells.tsv
    touch clonomap_out/families.tsv clonomap_out/reference_candidates.tsv clonomap_out/hc_unassigned.tsv
    """
}
