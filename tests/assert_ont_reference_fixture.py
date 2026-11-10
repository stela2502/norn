#!/usr/bin/env python3
"""Assert biological truth of Norn's published Nelrune transcript-mode outputs."""
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

out = Path(sys.argv[1])
tx = 'SYNTH_SPLICE4-001_skip_chrSyntheticMouse1:751-930'
gtf = out / 'new_transcripts.gtf'
assert gtf.is_file(), f'Missing novel transcript GTF: {gtf}'
exons = []
for line in gtf.read_text().splitlines():
    cols = line.split('\t')
    if len(cols) == 9 and f'transcript_id "{tx}"' in cols[8] and cols[2] == 'exon':
        exons.append((int(cols[3]), int(cols[4])))
assert sorted(exons) == [(201, 380), (1301, 1480), (1851, 2030)], f'Wrong mutant exons: {exons}'

snp = out / 'new_transcript_snp_evidence.tsv'
assert snp.is_file(), f'Missing SNP report: {snp}'
with snp.open() as handle:
    rows = [r for r in csv.DictReader(handle, delimiter='\t') if r['transcript_id'] == tx]
by_umi = defaultdict(dict)
cells = set()
for r in rows:
    cells.add(r['cell_id'])
    key = r['snp_id']
    assert key not in by_umi[r['umi']], f'Duplicate SNP/UMI evidence: {r}'
    by_umi[r['umi']][key] = (int(r['ref_reads']), int(r['alt_reads']))
assert len(cells) == 1 and len(by_umi) == 2, f'Expected one cell and two mutant UMIs, got {cells}, {by_umi}'
for umi, calls in by_umi.items():
    assert calls.get('1') == (0, 1), f'{umi}: exon-3 reporter should be ALT: {calls}'
    assert calls.get('2') == (1, 0), f'{umi}: other reporter should be REF: {calls}'

regions = out / 'new_transcript_supported_regions.tsv'
assert regions.is_file(), f'Missing molecule regions: {regions}'
with regions.open() as handle:
    supported = {r['umi'] for r in csv.DictReader(handle, delimiter='\t') if r['transcript_id'] == tx}
assert supported == set(by_umi), f'Supported UMIs differ from SNP UMIs: {supported} vs {set(by_umi)}'

matrix_dir = out / 'new_transcripts_mtx'
with gzip.open(matrix_dir / 'features.tsv.gz', 'rt') as handle:
    features = [line.rstrip('\n').split('\t') for line in handle]
with gzip.open(matrix_dir / 'barcodes.tsv.gz', 'rt') as handle:
    barcodes = [line.strip() for line in handle]
feature_indices = [i + 1 for i, row in enumerate(features) if tx in row]
assert len(feature_indices) == 1, f'Novel transcript not uniquely present in matrix features: {feature_indices}'
assert len(barcodes) == 1, f'Expected one cell barcode in matrix: {barcodes}'
with gzip.open(matrix_dir / 'matrix.mtx.gz', 'rt') as handle:
    lines = (line.strip() for line in handle if line.strip() and not line.startswith('%'))
    dims = [int(x) for x in next(lines).split()]
    entries = [line.split() for line in lines]
assert dims[0] == len(features) and dims[1] == len(barcodes), f'Bad matrix dimensions: {dims}'
value = sum(float(v) for row, col, v in entries if int(row) == feature_indices[0] and int(col) == 1)
assert value == 2, f'Expected 2 mutant transcript UMIs in expression matrix, got {value}'
print(f'[norn fixture] PASS: mutant transcript {tx}, 2 UMIs, ALT/REF SNP evidence, expression=2')
