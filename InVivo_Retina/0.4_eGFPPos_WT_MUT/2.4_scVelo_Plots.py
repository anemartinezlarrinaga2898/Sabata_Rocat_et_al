#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate publication-oriented RNA-velocity stream plots for WT and MUT.

All historical plotting variants are retained for traceability. See
REVIEW_NOTES.md before deciding which variant(s) to keep in the final repository.
"""

import numpy as np
import scanpy as sc
import scvelo as scv
from pathlib import Path
import os

scv.settings.verbosity = 3

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

# Directory used by Scanpy/scVelo when save=... is specified
sc.settings.figdir = str(path_out)
# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
adata_wt = sc.read_h5ad(path_out / "AnnData_WT_VELOCITY.h5ad")
adata_mut = sc.read_h5ad(path_out / "AnnData_MUT_VELOCITY.h5ad")

cluster_order = [
    "Angiogenic EC",
    "Proliferative 1",
    "Proliferative 2",
    "VenousCapillary EC",
    "Capillary",
    "ArteryCapillary EC",
    "Arterial EC",
]

cluster_colors = {
    "Angiogenic EC": "#f4cedb",
    "Proliferative 1": "#b5dcfb",
    "Proliferative 2": "#90caf9",
    "VenousCapillary EC": "#C5E1A5",
    "Capillary": "#e1be8f",
    "ArteryCapillary EC": "#FFE082",
    "Arterial EC": "#FF8A65",
}

cond_name= "MUT"
scv.pl.velocity_embedding_stream(
    adata_mut,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    show=False
)

plt.savefig(
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify.pdf",
    format="pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# WT

cond_name= "WT"
scv.pl.velocity_embedding_stream(
    adata_wt,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    show=False
)

plt.savefig(
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V2.pdf",
    format="pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

from PIL import Image
import matplotlib.pyplot as plt

cond_name = "WT"

output_png = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V2.png"
)

output_pdf = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V2.pdf"
)

fig, ax = plt.subplots(figsize=(8, 7))

scv.pl.velocity_embedding_stream(
    adata_wt,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    size=300,
    ax=ax,
    show=False
)

legend = ax.get_legend()
if legend is not None:
    legend.remove()

fig.savefig(
    output_png,
    format="png",
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)

Image.open(output_png).convert("RGB").save(
    output_pdf,
    "PDF",
    resolution=600
)

print(f"PNG: {output_png}")
print(f"PDF: {output_pdf}")

# Historical plotting variant 

cond_name = "WT"

output_pdf = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V3.pdf"
)

# Replace non-finite embedding values before vector PDF export
if "X_umap" in adata_wt.obsm:
    adata_wt.obsm["X_umap"] = np.nan_to_num(adata_wt.obsm["X_umap"], nan=0.0, posinf=0.0, neginf=0.0)

if "velocity_umap" in adata_wt.obsm:
    adata_wt.obsm["velocity_umap"] = np.nan_to_num(adata_wt.obsm["velocity_umap"], nan=0.0, posinf=0.0, neginf=0.0)
# ----------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 7))

scv.pl.velocity_embedding_stream(
    adata_wt,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    size=300,          # Same point-size convention as MUT
    ax=ax,
    show=False
)

legend = ax.get_legend()
if legend is not None:
    legend.remove()

fig.savefig(
    output_pdf,
    format="pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)

print(f"PDF vectorial guardado sin errores en: {output_pdf}")

# Historical plotting variant 

from PIL import Image
import matplotlib.pyplot as plt

cond_name = "MUT"

output_png = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V10.png"
)

output_pdf = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V10.pdf"
)

# Figure dimensions
fig, ax = plt.subplots(figsize=(10, 10))

# Generate stream plot with matched parameters
scv.pl.velocity_embedding_stream(
    adata_mut,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    size=600,
    density=2.0,
    linewidth = 1.5,
    arrow_size=1.5,
    ax=ax,
    show=False
)

legend = ax.get_legend()
if legend is not None:
    legend.remove()

ax.set_title("")

# Save high-resolution PNG
fig.savefig(
    output_png,
    format="png",
    dpi=1000,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none"
)

plt.close(fig)

# Convert PNG to PDF
Image.open(output_png).convert("RGB").save(
    output_pdf,
    "PDF",
    resolution=1000
)

print(f"MUT PNG: {output_png}")
print(f"MUT PDF: {output_pdf}")

# -- WT 

from PIL import Image
import matplotlib.pyplot as plt

cond_name = "WT"

output_png = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V10.png"
)

output_pdf = (
    path_out
    / f"VelocityStream_UMAP_Layer_1_{cond_name}_Modify_V10.pdf"
)

# Figure dimensions
fig, ax = plt.subplots(figsize=(10, 10))

# Generate stream plot with matched parameters
scv.pl.velocity_embedding_stream(
    adata_wt,
    basis="umap",
    color="Layer_1",
    palette=cluster_colors,
    legend_loc="none",
    size=600,
    density=2.0,
    linewidth = 1.5,
    arrow_size=1.5,
    ax=ax,
    show=False
)

legend = ax.get_legend()
if legend is not None:
    legend.remove()

ax.set_title("")

# Save high-resolution PNG
fig.savefig(
    output_png,
    format="png",
    dpi=1000,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none"
)

plt.close(fig)

# Convert PNG to PDF
Image.open(output_png).convert("RGB").save(
    output_pdf,
    "PDF",
    resolution=1000
)

print(f"WT PNG: {output_png}")
print(f"WT PDF: {output_pdf}")
