# Norn — initial Lumrik orchestration prototype

Norn is a deliberately small Nextflow DSL2 workflow around the current Lumrik command-line programs. It is not yet an nf-core pipeline.

The normal Norn input is **raw reference material**, not a collection of manually prebuilt indexes:

```text
genome.fa + annotation.gtf
          |
          +--> BUILD_STAR_INDEX   --> results/reference/star/star_index/
          |
          +--> BUILD_SPLICE_INDEX --> results/reference/splice/reference.splice.idx
          |
          +--> BUILD_VDJ_INDEX    --> results/reference/vdj/reference.vdjidx
```

Those are ordinary Nextflow processes. They are therefore independently cached by `-resume`, and the generated reference assets are also published as first-class Norn results for reuse by later runs.

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
- `cell_barcode_len`: optional; use `27` for legacy BD MEX normalization in `nelrune-vdj`

## Normal first run

Only the genome and annotation need to be supplied as reference inputs:

```bash
nextflow run main.nf \
  -profile local \
  --samplesheet assets/samplesheet.example.csv \
  --gtf /refs/genes.gtf \
  --genome /refs/genome.fa \
  --mapper star \
  --mapper_threads 8 \
  --nelrune_threads 8 \
  --outdir results
```

The first run builds the STAR, Lumrik splice and VDJ references. Repeating the run with `-resume` reuses those tasks unless their actual inputs or process definitions changed:

```bash
nextflow run main.nf \
  -profile local \
  --samplesheet assets/samplesheet.example.csv \
  --gtf /refs/genes.gtf \
  --genome /refs/genome.fa \
  --outdir results \
  -resume
```

The published reference assets can also be reused explicitly in another Norn run with `--mapper_index`, `--splice_index`, and `--vdj_index`. These options are escape hatches, not the default workflow.

## First-class reference results

With the default `--outdir results`:

```text
results/
  reference/
    star/
      star_index/
    splice/
      reference.splice.idx
    vdj/
      reference.vdjidx
  <sample>/
    nelrune/
    vdj/
```

`reference_publish_mode='link'` and `publish_mode='link'` are the defaults, avoiding duplicate large BAM/index files when the Nextflow work directory and output directory share a filesystem. Set either mode to `copy` if hard links are unsuitable for the installation.

## Deliberate milestone-1 limits

- STAR is the only mapper for which Norn currently builds the mapper index itself. Other mapper backends can still be used by supplying `--mapper_index` until they get their own reference module.
- ONT/BAM sample mode is not wired yet.
- Additional-feature FASTA/built-in wiring and standalone `lumrik-guides` orchestration are not wired yet; current Nelrune still couples feature finalization with Beacon calling.
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
