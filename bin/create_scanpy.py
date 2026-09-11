#!/usr/bin/env python3

import argparse
import gzip
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a Scanpy/AnnData object from Norn outputs."
    )
    parser.add_argument("--exonic", required=True)
    parser.add_argument("--intronic", required=True)
    parser.add_argument("--project", default="Norn")
    parser.add_argument("--output", required=True)
    parser.add_argument("--vdj-calls")
    parser.add_argument("--vdj-receptors")
    return parser.parse_args()


def read_features(directory: str) -> pd.DataFrame:
    gz_path = Path(directory) / "features.tsv.gz"
    plain_path = Path(directory) / "features.tsv"
    path = gz_path if gz_path.exists() else plain_path
    if not path.exists():
        raise FileNotFoundError(f"Missing 10x feature table: {gz_path} or {plain_path}")

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        features = pd.read_csv(handle, sep="\t", header=None, dtype=str)

    if features.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")
    while features.shape[1] < 3:
        features[features.shape[1]] = pd.NA

    features = features.iloc[:, :3].copy()
    features.columns = ["feature_id", "feature_name", "feature_type"]
    features["gene_id"] = features["feature_id"]
    features["gene_name"] = features["feature_name"]
    return features


def make_unique(values):
    seen = {}
    result = []
    for raw in values:
        value = str(raw)
        n = seen.get(value, 0)
        result.append(value if n == 0 else f"{value}-{n}")
        seen[value] = n + 1
    return result


def read_10x(directory: str):
    # Lumrik emits gzipped 10x MEX files; the deliberately tiny nf-test
    # fixtures are plain text. Scanpy defaults to compressed Cell Ranger v3
    # input, so detect which representation is present.
    directory_path = Path(directory)
    compressed = (directory_path / "matrix.mtx.gz").exists()
    matrix_path = directory_path / ("matrix.mtx.gz" if compressed else "matrix.mtx")
    if not matrix_path.exists():
        raise FileNotFoundError(
            f"Missing 10x matrix: {directory_path / 'matrix.mtx.gz'} or "
            f"{directory_path / 'matrix.mtx'}"
        )

    # var_names='gene_symbols' mirrors the Seurat object: feature names are the
    # user-facing gene symbols while the stable IDs are retained in metadata.
    adata = sc.read_10x_mtx(
        directory,
        var_names="gene_symbols",
        make_unique=True,
        gex_only=True,
        compressed=compressed,
    )
    return adata


def align_matrix(source, source_features, source_cells, features, cells):
    source = sparse.coo_matrix(source)
    feature_index = pd.Index(features)
    cell_index = pd.Index(cells)
    fi = feature_index.get_indexer(pd.Index(source_features))
    ci = cell_index.get_indexer(pd.Index(source_cells))

    row_map = fi[source.row]
    col_map = ci[source.col]
    keep = (row_map >= 0) & (col_map >= 0)

    return sparse.csr_matrix(
        (source.data[keep], (row_map[keep], col_map[keep])),
        shape=(len(features), len(cells)),
    )


def combine_feature_metadata(primary, secondary, features):
    columns = list(dict.fromkeys(list(primary.columns) + list(secondary.columns)))
    out = pd.DataFrame(index=pd.Index(features, name=None), columns=columns, dtype="object")

    for src in (primary, secondary):
        common_rows = out.index.intersection(src.index)
        for column in src.columns:
            current = out.loc[common_rows, column]
            missing = current.isna() | current.astype(str).eq("")
            if missing.any():
                fill_rows = common_rows[missing.to_numpy()]
                out.loc[fill_rows, column] = src.loc[fill_rows, column].astype("object")
    return out


def productivity_flag(status):
    if pd.isna(status) or status == "" or status == "unknown_no_v_cds":
        return pd.NA
    return status == "productive"


def add_vdj_metadata(adata, calls_path, receptors_path):
    calls = pd.read_csv(calls_path, sep="\t", dtype=str)
    receptors = pd.read_csv(receptors_path, sep="\t", dtype=str)

    needed_calls = {"cell", "recombination_id", "productivity_status"}
    needed_receptors = {"cell", "heavy_recombination_id", "light_recombination_id"}
    missing_calls = needed_calls.difference(calls.columns)
    missing_receptors = needed_receptors.difference(receptors.columns)
    if missing_calls:
        raise ValueError(f"VDJ calls table lacks required columns: {', '.join(sorted(missing_calls))}")
    if missing_receptors:
        raise ValueError(
            f"VDJ receptor table lacks required columns: {', '.join(sorted(missing_receptors))}"
        )

    receptors = receptors.drop_duplicates("cell", keep="first").set_index("cell")
    cells = pd.Index(adata.obs_names)
    aligned = receptors.reindex(cells)

    hc = aligned["heavy_recombination_id"].replace("", pd.NA)
    lc = aligned["light_recombination_id"].replace("", pd.NA)

    lookup = {
        (row.cell, row.recombination_id): row.productivity_status
        for row in calls.itertuples(index=False)
    }

    hc_status = [lookup.get((cell, rid), pd.NA) if pd.notna(rid) else pd.NA for cell, rid in zip(cells, hc)]
    lc_status = [lookup.get((cell, rid), pd.NA) if pd.notna(rid) else pd.NA for cell, rid in zip(cells, lc)]

    adata.obs["HC_recombination_id"] = hc.array
    adata.obs["LC_recombination_id"] = lc.array
    adata.obs["HC_productive"] = pd.array([productivity_flag(x) for x in hc_status], dtype="boolean")
    adata.obs["LC_productive"] = pd.array([productivity_flag(x) for x in lc_status], dtype="boolean")


def main():
    args = parse_args()
    if bool(args.vdj_calls) != bool(args.vdj_receptors):
        raise ValueError("--vdj-calls and --vdj-receptors must be supplied together")

    print(f"Reading exonic matrix: {args.exonic}")
    exonic_ad = read_10x(args.exonic)
    exonic_meta = read_features(args.exonic)
    exonic_meta.index = exonic_ad.var_names

    print(f"Reading intronic matrix: {args.intronic}")
    intronic_ad = read_10x(args.intronic)
    intronic_meta = read_features(args.intronic)
    intronic_meta.index = intronic_ad.var_names

    # Match the Seurat object semantics: exonic cells are canonical; features
    # are the union of exon and intron features.
    cells = list(exonic_ad.obs_names)
    features = list(dict.fromkeys(list(exonic_ad.var_names) + list(intronic_ad.var_names)))

    exonic = align_matrix(
        exonic_ad.X.T,
        exonic_ad.var_names,
        exonic_ad.obs_names,
        features,
        cells,
    )
    intronic = align_matrix(
        intronic_ad.X.T,
        intronic_ad.var_names,
        intronic_ad.obs_names,
        features,
        cells,
    )

    var = combine_feature_metadata(exonic_meta, intronic_meta, features)
    obs = pd.DataFrame(index=pd.Index(cells, name=None))
    obs["cell_id"] = cells

    # AnnData is cells x genes, hence the transpose from Lumrik's 10x matrices.
    adata = ad.AnnData(X=exonic.T.tocsr(), obs=obs, var=var)
    adata.layers["counts"] = exonic.T.tocsr()
    adata.layers["spliced"] = exonic.T.tocsr()
    adata.layers["unspliced"] = intronic.T.tocsr()

    if args.vdj_calls:
        print("Adding VDJ cell metadata")
        add_vdj_metadata(adata, args.vdj_calls, args.vdj_receptors)
    else:
        adata.obs["HC_recombination_id"] = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
        adata.obs["LC_recombination_id"] = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
        adata.obs["HC_productive"] = pd.Series(pd.NA, index=adata.obs_names, dtype="boolean")
        adata.obs["LC_productive"] = pd.Series(pd.NA, index=adata.obs_names, dtype="boolean")

    adata.uns["norn"] = {
        "created_by": "Norn",
        "project": args.project,
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "RNA_layers": ["counts", "spliced", "unspliced"],
    }

    print(f"Writing AnnData object: {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Created {args.output}")
    print(f"  cells: {adata.n_obs}")
    print(f"  RNA features: {adata.n_vars}")
    named_layers = [key for key in adata.layers.keys() if key is not None]
    print(f"  layers: {', '.join(named_layers)}")


if __name__ == "__main__":
    main()
