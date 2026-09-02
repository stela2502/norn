process NELRUNE_VDJ {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/vdj", mode: params.publish_mode

    input:
    tuple val(meta), path(exonic), path(bam)
    path vdj_index

    output:
    tuple val(meta), path('vdj_out/vdj_calls.tsv'), path('vdj_out/vdj_clones.tsv'), path('vdj_out/vdj-mapping-info.txt'), emit: tables
    tuple val(meta), path('vdj_out/vdj_observed.fasta'), path('vdj_out/vdj_naive.fasta'), optional: true, emit: sequences

    script:
    def barcodeArg = meta.cell_barcode_len ? "--cell-barcode-len ${meta.cell_barcode_len}" : ''
    def shardArg = params.vdj_bam_shards ? "--bam-shards ${params.vdj_bam_shards}" : ''
    def seqArg = params.vdj_write_sequences ? '--write-sequences' : ''
    """
    nelrune-vdj \\
        --exonic ${exonic} \\
        --bam ${bam} \\
        --index ${vdj_index} \\
        --out vdj_out \\
        ${barcodeArg} \\
        ${shardArg} \\
        ${seqArg}
    """
}
