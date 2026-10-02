#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Transfer spliced/unspliced layers to the cleaned eGFP-positive AnnData objects.

The script matches barcodes and genes between the cleaned AnnData objects and
historical scVelo source objects, transfers splicing layers, and concatenates
selected WT and MUT samples.
"""

import os
from pathlib import Path

import anndata as ad
import scanpy as sc

PROJECT_ROOT = Path(os.environ.get("PROJECT_ROOT", ".")).resolve()
ANN_DATA_DIR = PROJECT_ROOT / "0.6_Join_MUT_WT" / "EC_Subset" / "eGFP" / "scVelo" / "AnnDataObj"
SCVELO_SOURCE_DIR = PROJECT_ROOT / "0.3_Dowstream_Analysis" / "TrayectoryAnalysis" / "AnnDataObj"
ANN_DATA_DIR.mkdir(parents=True, exist_ok=True)

path_guardar_tray = ANN_DATA_DIR
path_clean_data = ANN_DATA_DIR
path_scvelo_data = SCVELO_SOURCE_DIR

samples_to_use = ["L1153", "L1154", "L1152", "L1155"]  # Selected WT and MUT samples

print(f"\nTarget samples to process: {samples_to_use}")
###################################################################################################

# ------------------------------------------------------------------
# Process and match each sample
# ------------------------------------------------------------------

adata_list = []
patient_ids = []

for patient_id in samples_to_use:
    clean_path = path_clean_data / f"AnnData_{patient_id}.h5ad"
    scvelo_path = path_scvelo_data / f"AnnData_{patient_id}_scVelo.h5ad"
    
    if not clean_path.exists():
        print(f"⚠️ [{patient_id}] Clean file not found at {clean_path} -> skipping.")
        continue
    if not scvelo_path.exists():
        print(f"⚠️ [{patient_id}] scVelo file not found at {scvelo_path} -> skipping.")
        continue

    print(f"\n[{patient_id}] Loading clean AnnData and scVelo source...")

    # 1. Cargar ambos objetos
    adata_clean = sc.read_h5ad(clean_path)
    adata_scv = sc.read_h5ad(scvelo_path)

    # Harmonize cell barcodes
    # Remove the suffix after the first hyphen to match the historical scVelo objects
    adata_clean.obs_names = [bc.split('-')[0] for bc in adata_clean.obs_names]

    # Ensure unique cell and gene names
    adata_clean.obs_names_make_unique()
    adata_clean.var_names_make_unique()
    adata_scv.obs_names_make_unique()
    adata_scv.var_names_make_unique()

    # Identify cells and genes shared by both objects
    common_barcodes = adata_clean.obs_names.intersection(adata_scv.obs_names)
    common_genes = adata_clean.var_names.intersection(adata_scv.var_names)

    print(f"[{patient_id}] Matching: {len(common_barcodes)} cells and {len(common_genes)} genes in common.")

    # Match the clean object to the cells and genes available in the scVelo source object
    adata = adata_clean[common_barcodes, common_genes].copy()
    adata_scv = adata_scv[common_barcodes, common_genes].copy()

    # Transfer splicing layers from the historical scVelo source object
    for layer in ["spliced", "unspliced", "ambiguous"]:
        if layer in adata_scv.layers:
            adata.layers[layer] = adata_scv.layers[layer]

    # Store sample identifier
    adata.obs["patient"] = patient_id

    # Save sample-level object ready for velocity analysis
    merged_path = path_guardar_tray / f"AnnData_{patient_id}_ReadyForVelocity.h5ad"
    adata.write(merged_path)
    print(f"[{patient_id}] Successfully merged layers into clean annotations and saved.")

    adata_list.append(adata)
    patient_ids.append(patient_id)

# ------------------------------------------------------------------
# Combine selected patients into one final AnnData object
# ------------------------------------------------------------------

if len(adata_list) > 1:
    print("\nMerging selected patients into one combined AnnData...")

    ad_dict = {pid: ad_obj for pid, ad_obj in zip(patient_ids, adata_list)}

    adata_combined = ad.concat(
        ad_dict,
        join="inner",
        index_unique="-"
    )

    merged_all_path = path_guardar_tray / "AnnData_Combined_scVelo.h5ad"
    adata_combined.write(merged_all_path)
    print(f"Combined AnnData successfully saved → {merged_all_path}")

elif len(adata_list) == 1:
    print("\nWarning: only one valid sample was loaded.")
    adata_combined = adata_list[0]
    merged_all_path = path_guardar_tray / "AnnData_Combined_scVelo.h5ad"
    adata_combined.write(merged_all_path)

else:
    print("\nError: no samples could be processed.")

print("\nScript completed successfully.")



