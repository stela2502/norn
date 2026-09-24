process SC_ANALYSE {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_mode

    input:
    tuple val(meta), path(exonic)

    output:
    tuple val(meta), path('sc_analysis'), emit: results

    script:
    def excludeVdjArg = params.sc_analyse_exclude_vdj ? '--exclude-vdj-before-normalization' : ''
    """
    sc-analyse ${exonic} \
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
