include { BUILD_STAR_INDEX }   from './modules/local/build_star_index'
include { BUILD_SPLICE_INDEX } from './modules/local/build_splice_index'
include { BUILD_VDJ_INDEX }    from './modules/local/build_vdj_index'
include { NELRUNE_PREPARE }    from './modules/local/nelrune_prepare'
include { STAR_ALIGN }         from './modules/local/star'
include { NELRUNE_QUANT }      from './modules/local/nelrune_quant'
include { NELRUNE_VDJ }        from './modules/local/nelrune_vdj'

nextflow.enable.dsl=2

def requireFile(value, label) {
    if (!value) error("Missing required parameter: ${label}")
    file(value, checkIfExists: true)
}


def starIndexState(value) {
    if (!value) return [path: null, exists: false, empty: true, valid: false]

    def p = java.nio.file.Paths.get(value.toString()).toAbsolutePath().normalize()
    def exists = java.nio.file.Files.exists(p)
    def empty = true

    if (exists) {
        if (!java.nio.file.Files.isDirectory(p)) {
            return [path: p, exists: true, empty: false, valid: false]
        }
        def stream = java.nio.file.Files.list(p)
        try {
            empty = !stream.findAny().isPresent()
        } finally {
            stream.close()
        }
    }

    def required = ['Genome', 'SA', 'SAindex', 'genomeParameters.txt']
    def valid = exists && required.every { java.nio.file.Files.isRegularFile(p.resolve(it)) }
    [path: p, exists: exists, empty: empty, valid: valid]
}


def lumrikIndexState(value, runtimeId) {
    if (!value) return [path: null, exists: false, empty: true, valid: false, runtimeMatch: false]

    def p = java.nio.file.Paths.get(value.toString()).toAbsolutePath().normalize()
    def exists = java.nio.file.Files.exists(p)
    def regular = exists && java.nio.file.Files.isRegularFile(p)
    def size = regular ? java.nio.file.Files.size(p) : 0L
    def marker = java.nio.file.Paths.get(p.toString() + '.lumrik-runtime')
    def markerExists = java.nio.file.Files.isRegularFile(marker)
    def recordedRuntime = markerExists ? java.nio.file.Files.readString(marker).trim() : null
    def runtimeMatch = markerExists && recordedRuntime == runtimeId.toString()

    [path: p, marker: marker, exists: exists, empty: !exists || (regular && size == 0L),
     valid: regular && size > 0L && runtimeMatch, runtimeMatch: runtimeMatch,
     recordedRuntime: recordedRuntime]
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

            def meta = [id: id, chemistry: chemistry]
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

    // Mapper reference. For STAR, --mapper_index is a persistent index location:
    // reuse a complete index there, or build/publish one there when the path is
    // missing or empty. Without --mapper_index, keep the normal result location.
    if (params.mapper != 'star') {
        if (!params.mapper_index) {
            error("Norn can currently build mapper indexes only for --mapper star; supply --mapper_index for '${params.mapper}'")
        }
        mapper_idx_ch = Channel.value(requireFile(params.mapper_index, '--mapper_index'))
    } else {
        def requestedMapperIndex = params.mapper_index ?: "${params.outdir}/reference/star/star_index"
        def mapperState = starIndexState(requestedMapperIndex)

        if (mapperState.valid) {
            log.info "  STAR index: reusing ${mapperState.path}"
            mapper_idx_ch = Channel.value(file(mapperState.path.toString(), checkIfExists: true))
        } else {
            if (mapperState.exists && !mapperState.empty) {
                error("STAR index path exists but is incomplete: ${mapperState.path}. Remove/rename it or provide a valid STAR index.")
            }

            def indexPath = mapperState.path
            def indexParent = indexPath.parent
            if (indexParent == null) {
                error("Unable to determine parent directory for STAR index: ${indexPath}")
            }

            log.info "  STAR index: building ${indexPath}"
            BUILD_STAR_INDEX(
                genome_ch,
                gtf_ch,
                Channel.value(indexParent.toString()),
                Channel.value(indexPath.fileName.toString())
            )
            mapper_idx_ch = BUILD_STAR_INDEX.out.index
        }
    }

    // Lumrik indexes are persistent reference infrastructure, never analysis outputs.
    // Their configured locations are authoritative: reuse them when valid or rebuild
    // them in place when missing/stale, but never fall back to params.outdir.
    if (!params.splice_index) {
        error("--splice_index is required; provide the authoritative persistent Lumrik splice-index path")
    }
    def requestedSpliceIndex = params.splice_index
    def spliceState = lumrikIndexState(requestedSpliceIndex, params.container)
    if (spliceState.valid) {
        log.info "  splice index: reusing ${spliceState.path} (${spliceState.recordedRuntime})"
        splice_idx_ch = Channel.value(file(spliceState.path.toString(), checkIfExists: true))
    } else {
        if (spliceState.exists && !spliceState.empty && !spliceState.runtimeMatch) {
            log.info "  splice index: rebuilding ${spliceState.path} (missing/stale Lumrik runtime marker; current ${params.container})"
        } else if (spliceState.exists && !spliceState.empty) {
            error("Splice index path exists but is not a regular non-empty file: ${spliceState.path}. Remove/rename it or provide a valid splice index.")
        }
        def indexPath = spliceState.path
        def indexParent = indexPath.parent
        if (indexParent == null) error("Unable to determine parent directory for splice index: ${indexPath}")

        log.info "  splice index: building ${indexPath}"
        BUILD_SPLICE_INDEX(
            gtf_ch,
            Channel.value(indexParent.toString()),
            Channel.value(indexPath.fileName.toString()),
            Channel.value(params.container.toString())
        )
        splice_idx_ch = BUILD_SPLICE_INDEX.out.index
    }

    // VDJ reference follows the same authoritative persistent-index semantics.
    if (params.run_vdj) {
        if (!params.vdj_index) {
            error("--vdj_index is required when --run_vdj is enabled; provide the authoritative persistent Lumrik VDJ-index path")
        }
        def requestedVdjIndex = params.vdj_index
        def vdjState = lumrikIndexState(requestedVdjIndex, params.container)
        if (vdjState.valid) {
            log.info "  VDJ index: reusing ${vdjState.path} (${vdjState.recordedRuntime})"
            vdj_idx_ch = Channel.value(file(vdjState.path.toString(), checkIfExists: true))
        } else {
            if (vdjState.exists && !vdjState.empty && !vdjState.runtimeMatch) {
                log.info "  VDJ index: rebuilding ${vdjState.path} (missing/stale Lumrik runtime marker; current ${params.container})"
            } else if (vdjState.exists && !vdjState.empty) {
                error("VDJ index path exists but is not a regular non-empty file: ${vdjState.path}. Remove/rename it or provide a valid VDJ index.")
            }
            def indexPath = vdjState.path
            def indexParent = indexPath.parent
            if (indexParent == null) error("Unable to determine parent directory for VDJ index: ${indexPath}")

            log.info "  VDJ index: building ${indexPath}"
            BUILD_VDJ_INDEX(
                gtf_ch,
                genome_ch,
                Channel.value(indexParent.toString()),
                Channel.value(indexPath.fileName.toString()),
                Channel.value(params.container.toString())
            )
            vdj_idx_ch = BUILD_VDJ_INDEX.out.index
        }
    }

    NELRUNE_PREPARE(samples_ch)
    STAR_ALIGN(NELRUNE_PREPARE.out.prepared, mapper_idx_ch)

    quant_input_ch = STAR_ALIGN.out.bam
        .join(NELRUNE_PREPARE.out.prepared, by: 0)
        .map { meta, bam, prepare_dir -> tuple(meta, bam, prepare_dir) }

    NELRUNE_QUANT(quant_input_ch, splice_idx_ch)

    if (params.run_vdj) {
        vdj_input_ch = NELRUNE_QUANT.out.exonic
            .join(NELRUNE_QUANT.out.bam, by: 0)
            .map { meta, exonic, bam -> tuple(meta, exonic, bam) }

        NELRUNE_VDJ(vdj_input_ch, vdj_idx_ch)
    }


}
