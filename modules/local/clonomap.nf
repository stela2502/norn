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
    mkdir -p vdj_input
    ln -s "${vdj_calls}" vdj_input/vdj_calls.tsv
    ln -s "${airr}" vdj_input/airr_rearrangements.tsv
    clonomap_family --vdj-out vdj_input --out clonomap_out ${params.clonomap_extra_args ?: ''}
    """

    stub:
    """
    mkdir -p clonomap_out
    printf 'source\tcell\tfamily\tlc_clone\n' > clonomap_out/cells.tsv
    touch clonomap_out/families.tsv clonomap_out/reference_candidates.tsv clonomap_out/hc_unassigned.tsv
    """
}
