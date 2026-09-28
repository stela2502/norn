#!/usr/bin/env bash
set -euo pipefail

NORN="$HOME/git_Projects/norn/main.nf"
SAMPLESHEET="$HOME/NAS/NELRUNE/REAL_TEST/sample_sheet.csv"
OUTDIR="$HOME/NAS/NELRUNE/REAL_TEST_v0.7.0"

REFERENCE="$HOME/sens05_shared/common/genome/genomes/mouse/GRCm39.MM39"
GENOME="$REFERENCE/GRCm39.genome.fa"
GTF="$REFERENCE/gencode.vM39.chr_patch_hapl_scaff.annotation.gtf"

# Keep using the existing local indices.
# Do NOT rebuild/replace the indices used by the running server job.
INDEX_ROOT="$HOME/NAS/NELRUNE/indexes/GRCm39_M39"
STAR_INDEX="$INDEX_ROOT/star"
SPLICE_INDEX="$INDEX_ROOT/mouse_GRCm39_M39.splice.idx"
VDJ_INDEX="$INDEX_ROOT/mouse_GRCm39_M39.vdjidx"

RESUME=()
# RESUME=(-resume)

DEBUG=(--max_reads 1000000)

nextflow run "$NORN" \
    "${RESUME[@]}" \
    "${DEBUG[@]}" \
    -config local-memory.config \
    -profile local,singularity \
    --samplesheet "$SAMPLESHEET" \
    --additional_features bd_sample_mouse \
    --genome "$GENOME" \
    --gtf "$GTF" \
    --mapper_index "$STAR_INDEX" \
    --splice_index "$SPLICE_INDEX" \
    --vdj_index "$VDJ_INDEX" \
    --mapper_threads 4 \
    --prepare_threads 20 \
    --nelrune_threads 20 \
    --vdj_threads 8 \
    --min_umi_counts 50 \
    --run_vdj true \
    --health_server true \
    --nelrune_health_port_base 8787 \
    --vdj_health_port_base 18787 \
    --make_seurat false \
    --make_scanpy false \
    --outdir "$OUTDIR"
