#!/usr/bin/env bash
# Real, non-stub Norn ONT test against the Lumrik reference-backed experiment.
# Requires: nextflow, samtools, python3, and a Nelrune binary with materialize-norn-test-data.
# For containers, set NORN_PROFILE=local,apptainer or local,docker.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE="$ROOT/tests/data/ont_reference_fixture"
OUT="${NORN_FIXTURE_OUT:-$ROOT/target/ont-reference-fixture}"
NELRUNE_BIN="${NELRUNE_BIN:-nelrune}"
NORN_PROFILE="${NORN_PROFILE:-local,apptainer}"
: "${LUMRIK_IMAGE:?Set LUMRIK_IMAGE to the newly published ghcr.io/stela2502/lumrik:v... image}"
mkdir -p "$OUT" "$OUT/generated" "$OUT/indexes" "$OUT/results" "$OUT/fixture_source"
cp "$FIXTURE/experiment.yaml" "$OUT/fixture_source/experiment.yaml"
for name in genome.fa.gz genes.gtf.gz variants.vcf.gz; do
    base64 -d "$FIXTURE/$name.b64" > "$OUT/fixture_source/$name"
done

"$NELRUNE_BIN" materialize-norn-test-data --config "$OUT/fixture_source/experiment.yaml" --out "$OUT/generated"

# Preserve the actual 43 sequences and read IDs; create a valid unaligned BAM.
python3 - "$OUT/generated/ont_reads.fastq" "$OUT/reads.sam" <<'PYTHON'
import sys
fastq, sam = sys.argv[1:]
with open(fastq) as src, open(sam, 'w') as dst:
    dst.write('@HD\tVN:1.6\tSO:unknown\n')
    count = 0
    while True:
        header = src.readline()
        if not header:
            break
        seq, plus, qual = src.readline(), src.readline(), src.readline()
        assert header.startswith('@') and plus.startswith('+') and len(seq.strip()) == len(qual.strip())
        dst.write(f'{header[1:].strip()}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq.strip()}\t{qual.strip()}\n')
        count += 1
    assert count == 43, f'Expected 43 synthetic reads, got {count}'
PYTHON
samtools view -b -o "$OUT/ont_reads.bam" "$OUT/reads.sam"
[[ "$(samtools view -c "$OUT/ont_reads.bam")" == 43 ]]
printf 'sample,bam,chemistry\nont_fixture,%s,tenx-three-prime-v3\n' "$OUT/ont_reads.bam" > "$OUT/samplesheet.csv"

cd "$ROOT"
# Exercise both ONT mapper branches with the same BAM and biological truth.
# Each mapper has its own output and index directory; the Lumrik splice index is shared.
for mapper in minimap2 star; do
    echo "[norn fixture] Running ONT mapper: $mapper"
    mapper_index="$OUT/indexes/genome.splice.mmi"
    if [[ "$mapper" == star ]]; then
        mapper_index="$OUT/indexes/star"
    fi
    nextflow run main.nf -profile "$NORN_PROFILE" -c "$ROOT/tests/ont_reference_resources.config" \
        --container "$LUMRIK_IMAGE" \
        --samplesheet "$OUT/samplesheet.csv" \
        --ont true \
        --mapper "$mapper" \
        --genome "$OUT/generated/genome.fa.gz" \
        --gtf "$OUT/generated/genes.gtf.gz" \
        --vcf "$OUT/generated/variants.vcf.gz" \
        --splice_index "$OUT/indexes/reference.splice.idx" \
        --mapper_index "$mapper_index" \
        --star_index_extra_args '--genomeSAindexNbases 3 --genomeChrBinNbits 12' \
        --primer_structure 'TYPE:GEX+FIXED:CTACACGACGCTCTTCCGATCT:mm=2+TENX_CELL:3p-v3+UMI:12' \
        --quant_mode transcript \
        --min_umi_counts 1 \
        --run_vdj false \
        --skip_sc_analyse true \
        --health_server false \
        --mapper_threads 2 \
        --nelrune_quant_threads 2 \
        --outdir "$OUT/results/$mapper" "$@"

    python3 "$ROOT/tests/assert_ont_reference_fixture.py" \
        "$OUT/results/$mapper/ont_fixture/nelrune/nelrune_out"
done
printf '[norn fixture] PASS: minimap2 and STAR recovered the expected transcript, SNPs and UMIs\n'
