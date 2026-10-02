#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimate dynamical RNA velocity separately in WT and MUT eGFP-positive cells.

The historical analytical parameters are preserved, including filtering,
normalization, dynamic model fitting, velocity graph construction, latent time,
velocity pseudotime and CellRank VelocityKernel computation.
"""

import numpy as np
import scanpy as sc
import scvelo as scv
import cellrank as cr
from pathlib import Path
import os

scv.settings.verbosity = 3
cr.settings.verbosity = 2

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.dpi"] = 300
# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path_guardar_original = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()

path_out = path_guardar_original / "0.6_Join_MUT_WT" / "EC_Subset" / "eGFP" / "scVelo" /  "cellrank_vk_only"
path_out.mkdir(parents=True, exist_ok=True)

path_ann_data = path_guardar_original / "0.6_Join_MUT_WT" / "EC_Subset" / "eGFP" / "scVelo" / "AnnDataObj"

# Para que scvelo/scanpy guarden con `save=...`
sc.settings.figdir = str(path_out)
# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
adata = sc.read_h5ad(path_ann_data / "AnnData_Combined_scVelo.h5ad")

# ------------------------------------------------------------------
# QC quick plots (Globales)
# ------------------------------------------------------------------
# sc.pl.proportions(adata, groupby="Layer_1", save="_proportions_by_Layer_1.pdf")
# sc.pl.umap(adata, color="Layer_1", frameon=False, save="_UMAP_Layer_1.pdf")

# # Check for unspliced per cell 
# adata.obs["unspliced_counts"] = np.array(adata.layers["unspliced"].sum(axis=1)).flatten()
# sc.pl.violin(adata, keys="unspliced_counts", groupby="Layer_1", rotation=90, save="_Vln_Unspliced.pdf")

# ------------------------------------------------------------------
# Joint preprocessing (shared PCA/neighborhood space)
# ------------------------------------------------------------------
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)

# ------------------------------------------------------------------
# scVelo dynamics and velocity estimated separately for WT and MUT
# ------------------------------------------------------------------
genotype_col = 'Phenotype' 

# 1. Split objects by condition
adata_wt = adata[adata.obs[genotype_col] == 'WT'].copy()
adata_mut = adata[adata.obs[genotype_col] == 'MUT'].copy()

# 2. Back up Layer_1 annotations before scVelo preprocessing
annot_wt_backup = adata_wt.obs['Layer_1'].copy()
annot_mut_backup = adata_mut.obs['Layer_1'].copy()

# ------------------------------------------------------------------
# WT Analysis
# ------------------------------------------------------------------
print("--- Estimating dynamics and velocity for WT ---")
adata_wt
scv.pp.filter_and_normalize(adata_wt, min_shared_counts=10, n_top_genes=2000,enforce=True)

# Restore Layer_1 annotation and remove unused categories/colors
adata_wt.obs['Layer_1'] = annot_wt_backup.loc[adata_wt.obs_names]
adata_wt.obs['Layer_1'] = adata_wt.obs['Layer_1'].astype('category').cat.remove_unused_categories()
if 'Layer_1_colors' in adata_wt.uns:
    del adata_wt.uns['Layer_1_colors']

# scVelo pipeline
scv.tl.recover_dynamics(adata_wt, n_jobs=8)
scv.tl.velocity(adata_wt, mode="dynamical")
scv.tl.velocity_graph(adata_wt)
scv.tl.latent_time(adata_wt)
scv.tl.velocity_pseudotime(adata_wt)

# ------------------------------------------------------------------
# MUT Analysis
# ------------------------------------------------------------------
print("--- Estimating dynamics and velocity for MUT ---")
scv.pp.filter_and_normalize(adata_mut, min_shared_counts=10, n_top_genes=2000,enforce=True)

# Restore Layer_1 annotation and remove unused categories/colors
adata_mut.obs['Layer_1'] = annot_mut_backup.loc[adata_mut.obs_names]
adata_mut.obs['Layer_1'] = adata_mut.obs['Layer_1'].astype('category').cat.remove_unused_categories()
if 'Layer_1_colors' in adata_mut.uns:
    del adata_mut.uns['Layer_1_colors']

# scVelo pipeline
scv.tl.recover_dynamics(adata_mut, n_jobs=8)
scv.tl.velocity(adata_mut, mode="dynamical")
scv.tl.velocity_graph(adata_mut)
scv.tl.latent_time(adata_mut)
scv.tl.velocity_pseudotime(adata_mut)

# ------------------------------------------------------------------
# Save processed condition-specific objects
# ------------------------------------------------------------------
adata_wt.write_h5ad(path_out / "AnnData_WT_VELOCITY.h5ad")
adata_mut.write_h5ad(path_out / "AnnData_MUT_VELOCITY.h5ad")

# ------------------------------------------------------------------
# CellRank: VelocityKernel only (condition-specific analysis)
# ------------------------------------------------------------------
conditions = { "MUT": adata_mut, "WT": adata_wt}

for cond_name, adata_sub in conditions.items():
    print(f"\n==============================================")
    print(f"Running CellRank for: {cond_name}")
    print(f"==============================================")

    vk = cr.kernels.VelocityKernel(adata_sub)
    vk.compute_transition_matrix()  # OK

    vk.write_to_adata()
    adata_sub.write_h5ad(path_out / f"AnnData_with_kernel_{cond_name}.h5ad", compression="gzip")
    
    # ------------------------------------------------------------------
    # 1. scVelo plots: streamlines, latent time and velocity pseudotime
    # ------------------------------------------------------------------
    scv.pl.velocity_embedding_stream(
        adata_sub, basis="umap", color="Layer_1", legend_loc="right", 
        save=f"_VelocityStream_UMAP_Layer_1_{cond_name}.pdf"
    )
    
    scv.pl.scatter(
        adata_sub, color="latent_time", basis="umap", 
        save=f"_latent_time_umap_{cond_name}.pdf"
    )

    scv.pl.scatter(
        adata_sub, color='velocity_pseudotime', cmap='gnuplot', 
        save=f"_Pseudotime_Velocity_{cond_name}.pdf"
    )
    