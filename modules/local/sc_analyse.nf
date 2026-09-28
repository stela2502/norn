process SC_ANALYSE {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(exprs), path(clonomap_input)

    output:
    tuple val(meta), path('sc_analysis'), emit: results

    script:
    def excludeVdjArg = params.sc_analyse_exclude_vdj ? '--exclude-vdj-before-normalization' : ''
    """
    mkdir -p analysis_input/filtered
    ln -s "\$(readlink -f "${exprs}")" analysis_input/filtered/exprs
    if [[ -f "${clonomap_input}" ]]; then
        mkdir -p analysis_input/clonomap_out
        ln -s "\$(readlink -f "${clonomap_input}")" analysis_input/clonomap_out/cells.tsv
    fi
    sc-analyse analysis_input \
        --out sc_analysis \
        --min-umi-count ${params.min_umi_counts} \
        ${excludeVdjArg}
    """

    stub:
    """
    mkdir -p sc_analysis
    touch sc_analysis/.stub
    """
}
