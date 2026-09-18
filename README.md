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
    +--> exonic MEX --------+
    |                       |
    +--> mapper BAM --------+--> NELRUNE_VDJ <---- VDJ index
    |                             |
    |                             +--> vdj_calls.tsv
    |                             +--> vdj_receptors.tsv
    |                             +--> airr_rearrangements.tsv
    |
    +--> intronic MEX
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
  --mapper_threads 8 \
  --nelrune_threads 8 \
  --vdj_threads 8 \
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

For BD chemistries, Norn derives the matching `nelrune-vdj --bd-cell-version` automatically from the samplesheet chemistry (`bd-v1`, `bd-v2-96`, or `bd-v2-384`). This keeps the canonical corrected 27-base barcode and the positional BD/Rustody cell id in the same namespace without carrying the old `cell_barcode_len` workaround.

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

Run the real Seurat and Scanpy component tests with a container runtime, for example:

```bash
nf-test test --tag component --profile +singularity
```

Use `+docker` or `+apptainer` instead when appropriate. The component tests use tiny 3-cell 10x-style exon/intron matrices and VDJ tables under `tests/data/`.
