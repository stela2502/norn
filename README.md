# Norn — initial Lumrik orchestration prototype

Norn is a deliberately small Nextflow DSL2 workflow around the current Lumrik command-line programs. It is not yet an nf-core pipeline.

The normal Norn input is raw reference material plus authoritative persistent locations for the Lumrik indexes:

```text
genome.fa + annotation.gtf
          |
          +--> BUILD_STAR_INDEX   --> results/reference/star/star_index/
          |
          +--> BUILD_SPLICE_INDEX --> --splice_index
          |
          +--> BUILD_VDJ_INDEX    --> --vdj_index
```

The Lumrik splice and VDJ indexes are reference infrastructure, not analysis results. Norn reuses or rebuilds them at the explicitly configured paths and never creates fallback copies below `--outdir`.

## Implemented sample graph

```text
samples.csv
    |
 FASTQ lanes
    |
    v
 NELRUNE <--------- STAR index + Lumrik splice index
    |
    +--> filtered/exonic MEX +
    |                       |
    +--> mapper BAM --------+--> NELRUNE_VDJ <---- VDJ index
    |                             |
    |                             +--> vdj_calls.tsv
    |                             +--> vdj_receptors.tsv
    |                             +--> airr_rearrangements.tsv
    |
    +--> filtered/intronic MEX
    +--> raw/{exonic,intronic,...} MEX
    +--> Nelrune QC/log/metrics
```

## Samplesheet

CSV columns:

- `sample`: unique sample id
- `r1`: one or more R1 files separated by `;`
- `r2`: matching R2 files separated by `;`
- `chemistry`: exact Lumrik chemistry value

## Normal first run

Supply the genome and annotation plus authoritative persistent paths for the Lumrik splice and VDJ indexes:

```bash
nextflow run main.nf \
  -profile local \
  --samplesheet assets/samplesheet.example.csv \
  --gtf /refs/genes.gtf \
  --genome /refs/genome.fa \
  --splice_index /refs/lumrik/reference.splice.idx \
  --vdj_index /refs/lumrik/reference.vdjidx \
  --mapper star \
  --mapper_threads 18 \
  --nelrune_threads 18 \
  --vdj_threads 18 \
  --outdir results
```

The first run builds the STAR, Lumrik splice and VDJ references. Repeating the run with `-resume` reuses those tasks unless their actual inputs or process definitions changed:

```bash
nextflow run main.nf \
  -profile local \
  --samplesheet assets/samplesheet.example.csv \
  --gtf /refs/genes.gtf \
  --genome /refs/genome.fa \
  --splice_index /refs/lumrik/reference.splice.idx \
  --vdj_index /refs/lumrik/reference.vdjidx \
  --outdir results \
  -resume
```

`--splice_index` is required. `--vdj_index` is required when `--run_vdj` is enabled. These paths are the authoritative Lumrik indexes across runs; changing `--outdir` does not create or select another copy.

For BD chemistries, Norn derives the matching `nelrune vdj --bd-cell-version` automatically from the samplesheet chemistry (`bd-v1`, `bd-v2-96`, or `bd-v2-384`). This keeps the canonical corrected 27-base barcode and the positional BD/Rustody cell id in the same namespace without carrying the old `cell_barcode_len` workaround.

## Reference locations

STAR may still be built below `--outdir` when no `--mapper_index` is supplied. Lumrik splice and VDJ indexes are different: they live only at the configured `--splice_index` and `--vdj_index` paths. They are never published below the analysis results directory.

## Deliberate milestone-1 limits

- STAR is the only mapper for which Norn currently builds the mapper index itself. Other mapper backends can still be used by supplying `--mapper_index` until they get their own reference module.
- ONT/BAM sample mode is not wired yet.
- Additional-feature FASTA/built-in wiring and standalone `lumrik-guides` orchestration are not wired yet.
- No downstream R/Seurat/reporting yet.
- No nf-core schema/templates/modules yet.
- Container configuration hooks exist, but no canonical Lumrik container URI is invented here.

The important milestone-1 property is that **reference generation, expensive primary mapping, and VDJ analysis are distinct cache boundaries**.


## Containerized execution

Norn uses the Lumrik runtime image declared by `params.container`. The intended
production image contains STAR plus static Lumrik musl binaries. For Singularity
on an HPC system:

    nextflow run main.nf -profile slurm,singularity \
      --samplesheet samples.csv \
      --genome /path/reference.fa \
      --gtf /path/annotation.gtf \
      -resume

Pin `params.container` to a released image tag (or, for maximum reproducibility,
an OCI digest) rather than `latest`. Public GHCR images can be pulled anonymously
by Singularity/Nextflow.

## Testing

Norn uses [nf-test](https://www.nf-test.com/) with the same process/pipeline testing model used by nf-core.

Fast wiring smoke test (no bioinformatics tools or containers required):

```bash
nf-test test --tag stub
```


## Reference-backed ONT transcript integration test

The non-stub `tests/run_ont_reference_fixture.sh` reuses the 43-read synthetic
Lumrik experiment (`tests/data/ont_reference_fixture/`). It builds the
minimap2, STAR and Lumrik splice indexes, runs Norn with `--quant_mode transcript`
and a VCF, then checks the **published** Nelrune outputs: novel exon chain,
cell/UMI-resolved SNP evidence, supported molecules and expression count 2.
Both ONT mapper branches (minimap2 and STAR) are checked against the same
novel-transcript, SNP/UMI and matrix-expression truth. STAR indexing uses
small-genome parameters only in this test. Use `--ont true` with `--mapper star`
for ONT BAM inputs; existing paired-FASTQ STAR behavior is unchanged.
The test skips downstream `sc-analyse` because this one-cell fixture is for
transcript discovery, not clustering. Normal Norn runs still include it.

Requires `nextflow`, `samtools`, `python3`, and a local `nelrune` executable
for materializing the fixture. Use the newly released Lumrik image:

```bash
LUMRIK_IMAGE=ghcr.io/stela2502/lumrik:<new-tag> \
  NELRUNE_BIN=~/git_Projects/lumrik/target/release/nelrune \
  NORN_PROFILE=local,apptainer \
  bash tests/run_ont_reference_fixture.sh
```

Set `NORN_PROFILE=local,docker` when using Docker instead. Generated assets
are retained under `target/ont-reference-fixture/` and excluded from Git.
