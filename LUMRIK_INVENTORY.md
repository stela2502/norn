# Lumrik source inventory for Norn (2026-09-02 upload)

## Initial production stages

### nelrune
Primary single-cell processing. Illumina uses paired `--r1`/`--r2` lists; ONT uses `--bam` instead. Requires chemistry/primer configuration, a mapper and `--mapper-index`, and Lumrik splice `--index`. Optional genome/VCF quantification inputs are inherited from bam-tide. Output root contains `exonic/`, `intronic/`, optional `ref/` + `alt/`, `nelrune.mapper.bam` (unless overridden with `--bam-out`), `nelrune-report.txt`, `nelrune.log`, metrics TSV, normalization intermediates, and per-additional-feature directories. Additional-feature finalization currently runs sc-beacon internally.

### gtf-splice-index
Reference preparation. `gtf-splice-index build --annotation <GTF/GFF> --index <FILE>` creates the reusable splice index consumed by Nelrune/bam-tide. Also has stats/transcriptome utilities.

### vdj-index
Reference preparation. `vdj-index --gtf <GTF> --genome <FASTA> --out <FILE.vdjidx>` builds the preferred VDJ index.

### nelrune-vdj
Downstream VDJ production stage. Requires `--exonic <MEX_DIR>`, `--bam <BAM>`, and preferably `--index <VDJIDX>` (or GTF+genome directly). Writes `vdj_calls.tsv`, `vdj_clones.tsv`, `vdj-mapping-info.txt`, and optional observed/naive FASTA. `--cell-barcode-len 27` supports legacy BD normalization. BAM sharding is automatic (~128 cells/shard, clamped 1..64) or set via `--bam-shards`.

### lumrik-guides
Standalone downstream guide caller. Requires raw and filtered 10x-style feature matrices, output dir, optional feature type and threads. Writes ambient/model/call/assignment/statistics tables plus posterior MEX and log. In the current Nelrune path, comparable Beacon analysis is already invoked internally for each additional feature type.

## Diagnostics / utilities, not initial Norn stages

`vdj-summary`, `vdj-rich-cell`, and `vdj-decode` are diagnostic/forensic tools. `identify_primers`, `mapper-wrapper`, `bd-fast-map`, and most bam-tide binaries expose lower-level functionality already orchestrated by Nelrune. They are useful for debugging or specialized future workflows, but should not be inserted into the initial Norn graph merely because they are public binaries.

bam-tide also contains specialized utilities (`bam-coverage`, `bam-quant`, `bam-subset-tag`, `bam-transcriptome-to-genome`, ONT/Illumina normalizers, etc.). These are candidates only when Norn gains workflows that genuinely require them independently.

## Current file contracts

- GTF -> `gtf-splice-index` -> serialized splice index -> Nelrune `--index`.
- GTF + genome FASTA -> `vdj-index` -> `.vdjidx` -> `nelrune-vdj --index`.
- Nelrune `exonic/` is a standard MatrixMarket directory (`barcodes.tsv[.gz]`, `features.tsv[.gz]`/`genes.tsv[.gz]`, `matrix.mtx[.gz]`) and is read directly by `nelrune-vdj`.
- Nelrune retained mapper BAM carries cell/UMI identity understood by `NelruneIdentityResolver` and is read directly by `nelrune-vdj`.
- Nelrune additional-feature directories are 10x-style sparse matrices; Beacon results are nested below `<feature_type>/beacon/` in the integrated path.

## Workflow-boundary issues worth patching separately

1. Nelrune couples additional-feature matrix writing and Beacon calling in `FeatureTagCounts::finalize_and_write`. A switch such as `--no-feature-calling` (or a clean raw/filtered feature-output contract) would let Norn cache feature extraction independently from calling/model changes.
2. `nelrune-vdj` has no threads CLI. Its BAM sharding bounds memory, but Nextflow cannot currently express/use requested CPU parallelism meaningfully for the computation itself.
3. Nelrune accepts an optional `--bam-out` inherited from bam-tide, but a pipeline should keep the default output inside the task directory. Allowing arbitrary external output paths undermines Nextflow staging/caching, so Norn should not expose that option initially.
4. Several outputs are directory contracts rather than named manifest objects. This is workable in Nextflow, but a small machine-readable run manifest from Nelrune would make future schema/provenance validation easier.
5. The binary file is `gtf-splice-index`, while Clap advertises command name `splice-index`; harmless, but slightly confusing in help/logs/container smoke tests.
6. `nelrune-vdj` output currently contains structural IDs/codes and deletion/P/N measurements, but not yet all of the richer biological columns described for the intended direction (for example explicit germline/observed/naive sequences in the TSV). Norn should consume the stable tables rather than parse internal serialization.

## Norn reference-preparation decision

For milestone 1, Norn treats `genome.fa` and `annotation.gtf` as the normal reference inputs and builds/publishes three independent reusable assets:

- STAR genome index (`STAR --runMode genomeGenerate`)
- Lumrik splice index (`gtf-splice-index build`)
- Lumrik VDJ index (`vdj-index`)

Prebuilt index parameters remain optional overrides, but they are not the recommended/default route. This keeps reference construction reproducible and independently resumable while making the resulting assets reusable outside the run that created them.
