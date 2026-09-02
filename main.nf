include { BUILD_STAR_INDEX }   from './modules/local/build_star_index'
include { BUILD_SPLICE_INDEX } from './modules/local/build_splice_index'
include { BUILD_VDJ_INDEX }    from './modules/local/build_vdj_index'
include { NELRUNE }            from './modules/local/nelrune'
include { NELRUNE_VDJ }        from './modules/local/nelrune_vdj'

nextflow.enable.dsl=2

def requireFile(value, label) {
    if (!value) error("Missing required parameter: ${label}")
    file(value, checkIfExists: true)
}

def sampleChannel(samplesheet) {
    Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def id = row.sample?.trim()
            if (!id) error('samplesheet contains an empty sample value')

            def chemistry = row.chemistry?.trim()
            if (!chemistry) error("sample ${id}: chemistry is required")

            def r1s = row.r1.split(';').collect { file(it.trim(), checkIfExists: true) }
            def r2s = row.r2.split(';').collect { file(it.trim(), checkIfExists: true) }
            if (r1s.size() != r2s.size()) error("sample ${id}: r1/r2 lane counts differ")

            def cellLen = row.cell_barcode_len?.trim() ? row.cell_barcode_len.trim() as Integer : null
            def meta = [id: id, chemistry: chemistry, cell_barcode_len: cellLen]
            tuple(meta, r1s, r2s)
        }
}

workflow {
    if (!params.samplesheet) error('Use --samplesheet samples.csv')

    // Genome + annotation are Norn's primary reference inputs. By default Norn
    // turns them into all reusable indexes and publishes those as first-class results.
    def gtf = requireFile(params.gtf, '--gtf')
    def genome = requireFile(params.genome, '--genome')

    samples_ch = sampleChannel(params.samplesheet)
    gtf_ch = Channel.value(gtf)
    genome_ch = Channel.value(genome)

    // Mapper reference. STAR is the built-in/default reference builder for milestone 1.
    // An externally prepared mapper index remains an escape hatch for other mappers.
    if (params.mapper_index) {
        mapper_idx_ch = Channel.value(requireFile(params.mapper_index, '--mapper_index'))
    } else {
        if (params.mapper != 'star') {
            error("Norn can currently build mapper indexes only for --mapper star; supply --mapper_index for '${params.mapper}'")
        }
        BUILD_STAR_INDEX(genome_ch, gtf_ch)
        mapper_idx_ch = BUILD_STAR_INDEX.out.index
    }

    // Lumrik splice index.
    if (params.splice_index) {
        splice_idx_ch = Channel.value(requireFile(params.splice_index, '--splice_index'))
    } else {
        BUILD_SPLICE_INDEX(gtf_ch)
        splice_idx_ch = BUILD_SPLICE_INDEX.out.index
    }

    // VDJ reference is built independently of sample processing so it can be cached/reused.
    if (params.vdj_index) {
        vdj_idx_ch = Channel.value(requireFile(params.vdj_index, '--vdj_index'))
    } else if (params.run_vdj) {
        BUILD_VDJ_INDEX(gtf_ch, genome_ch)
        vdj_idx_ch = BUILD_VDJ_INDEX.out.index
    }

    NELRUNE(samples_ch, splice_idx_ch, mapper_idx_ch)

    if (params.run_vdj) {
        vdj_input_ch = NELRUNE.out.exonic
            .join(NELRUNE.out.bam, by: 0)
            .map { meta, exonic, bam -> tuple(meta, exonic, bam) }

        NELRUNE_VDJ(vdj_input_ch, vdj_idx_ch)
    }
}
