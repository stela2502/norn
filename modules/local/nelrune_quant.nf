process NELRUNE_QUANT {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}/nelrune", mode: params.publish_mode

    input:
    tuple val(meta), path(bam), path(prepare_dir)
    path splice_index
    path genome

    output:
    tuple val(meta), path('nelrune_out/exonic'), emit: exonic
    tuple val(meta), path('nelrune_out/intronic'), emit: intronic
    tuple val(meta), path('nelrune_out/nelrune-report.txt'), emit: qc

    script:
    def arg = { flag, value -> value != null && value.toString() != '' ? "${flag} ${value}" : '' }
    def boolArg = { flag, value -> value as boolean ? flag : '' }
    def listArg = { flag, value ->
        if (value == null || value.toString().trim() == '') return ''
        def vals = value instanceof Collection ? value : value.toString().split(',').collect { it.trim() }.findAll { it }
        vals ? "${flag} ${vals.join(' ')}" : ''
    }

    def primerStructureArg = arg('--primer-structure', params.primer_structure)
    def whitelistArg = arg('--whitelist', params.whitelist)
    // The genome FASTA is a SNP-only runtime dependency. Norn still stages it
    // through a path input so that, when SNP collection is requested, the
    // container sees the same reference that Norn used to build its indexes.
    def collectSnps = params.vcf != null && params.vcf.toString().trim() != ''
    def genomeArg = collectSnps ? "--genome ${genome}" : ''
    def vcfArg = arg('--vcf', params.vcf)
    def maxReadsArg = arg('--max-reads', params.max_reads)
    def readTagArg = listArg('--read-tag-table', params.read_tag_table)
    def bamOutArg = arg('--bam-out', params.bam_out)
    def featuresArg = listArg('--additional-features', params.additional_features)
    def minCellCountsArg = arg('--min-cell-counts', params.min_umi_counts)

    """
    mkdir -p lumrik_tmp
    export LUMRIK_TMPDIR="\$PWD/lumrik_tmp"

    nelrune quant \\
        --bam ${bam} \\
        --prepare ${prepare_dir} \\
        --chemistry ${meta.chemistry} \\
        ${primerStructureArg} \\
        ${whitelistArg} \\
        --whitelist-mismatches ${params.whitelist_mismatches} \\
        ${boolArg('--detect-reverse-complement', params.detect_reverse_complement)} \\
        --index ${splice_index} \\
        ${genomeArg} \\
        ${vcfArg} \\
        --quant-mode ${params.quant_mode} \\
        --min-mapq ${params.min_mapq} \\
        ${maxReadsArg} \\
        ${boolArg('--read1-only', params.read1_only)} \\
        ${boolArg('--no-genome-refine', params.no_genome_refine)} \\
        ${boolArg('--require-strand', params.require_strand)} \\
        ${boolArg('--require-exact-junction-chain', params.require_exact_junction_chain)} \\
        --max-5p-overhang-bp ${params.max_5p_overhang_bp} \\
        --max-3p-overhang-bp ${params.max_3p_overhang_bp} \\
        --allowed-intronic-gap-size ${params.allowed_intronic_gap_size} \\
        --snp-min-anchor ${params.snp_min_anchor} \\
        ${readTagArg} \\
        --rt-read-id-column ${params.rt_read_id_column} \\
        --rt-cell-column ${params.rt_cell_column} \\
        --rt-cell-qual-column ${params.rt_cell_qual_column} \\
        --rt-umi-column ${params.rt_umi_column} \\
        --rt-umi-qual-column ${params.rt_umi_qual_column} \\
        --rt-original-read-id-column ${params.rt_original_read_id_column} \\
        --analysis-type ${params.analysis_type} \\
        --cell-tag ${params.cell_tag} \\
        --umi-tag ${params.umi_tag} \\
        ${bamOutArg} \\
        ${featuresArg} \\
        --additional-feature-min-hits ${params.additional_feature_min_hits} \\
        ${minCellCountsArg} \\
        --outpath nelrune_out \\
        ${params.nelrune_quant_extra_args ?: ''}
    """

    stub:
    """
    mkdir -p nelrune_out/exonic nelrune_out/intronic
    touch nelrune_out/exonic/matrix.mtx.gz nelrune_out/exonic/features.tsv.gz nelrune_out/exonic/barcodes.tsv.gz
    touch nelrune_out/intronic/matrix.mtx.gz nelrune_out/intronic/features.tsv.gz nelrune_out/intronic/barcodes.tsv.gz
    touch nelrune_out/nelrune-report.txt
    """
}
