#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create per-sample AnnData objects from Seurat-exported matrices and metadata.

Inputs are expected in the Conversion_Seurat_To_Python directory and must include
counts, metadata, gene names and PCA coordinates for each sample.
"""

import os
from pathlib import Path

import anndata
import pandas as pd
from scipy import io

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
INPUT_DIR = PROJECT_ROOT / "0.6_Join_MUT_WT" / "EC_Subset" / "eGFP" / "scVelo" / "Conversion_Seurat_To_Python"
OUTPUT_DIR = PROJECT_ROOT / "0.6_Join_MUT_WT" / "EC_Subset" / "eGFP" / "scVelo" / "AnnDataObj"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def generate_anndata(sample_id: str, input_dir: Path, output_dir: Path):
    """Create and save an AnnData object for one sample."""
    counts_path = input_dir / f"counts_{sample_id}.mtx"
    x = io.mmread(counts_path)
    print(f"[{sample_id}] Matrix loaded: {x.shape}")

    metadata_csv = input_dir / f"metadata_{sample_id}.csv"
    metadata_xlsx = input_dir / f"metadata_{sample_id}.xlsx"
    if metadata_csv.exists():
        cell_meta = pd.read_csv(metadata_csv)
    elif metadata_xlsx.exists():
        cell_meta = pd.read_excel(metadata_xlsx)
    else:
        raise FileNotFoundError(f"No metadata file found for {sample_id}")
    print(f"[{sample_id}] Metadata loaded: {cell_meta.shape}")

    gene_path = input_dir / f"gene_names_{sample_id}.csv"
    gene_names = pd.read_csv(gene_path, header=None).squeeze().tolist()
    print(f"[{sample_id}] Gene names loaded: {len(gene_names)} genes")

    adata = anndata.AnnData(X=x.transpose().tocsr())
    adata.obs = cell_meta
    adata.obs.index = adata.obs["barcode"]
    adata.var.index = gene_names

    pca_path = input_dir / f"pca_{sample_id}.csv"
    pca = pd.read_csv(pca_path)
    pca.index = adata.obs.index
    adata.obsm["X_pca"] = pca.to_numpy()

    if {"UMAP_1", "UMAP_2"}.issubset(adata.obs.columns):
        adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy()
    else:
        print(f"[{sample_id}] Warning: UMAP_1/UMAP_2 columns not found in metadata.")

    output_path = output_dir / f"AnnData_{sample_id}.h5ad"
    adata.write(output_path)
    print(f"[{sample_id}] AnnData saved -> {output_path}")
    return adata


count_files = sorted(INPUT_DIR.glob("counts_*.mtx"))
for count_file in count_files:
    sample_id = count_file.stem.replace("counts_", "")
    print(f"\n=== Processing {sample_id} ===")
    try:
        generate_anndata(sample_id, INPUT_DIR, OUTPUT_DIR)
    except Exception as exc:
        print(f"[{sample_id}] ERROR: {exc}")
