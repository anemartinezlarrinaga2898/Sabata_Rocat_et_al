# Review notes — eGFP-positive WT + MUT analysis

These scripts were homogenized for repository readability without intentionally changing the scientific analysis. The following points should be reviewed before treating the folder as a fully reproducible final pipeline.

1. **`3.0_Subset_eGFP.R` contains a historical object-assignment inconsistency.** `SeuratPipeline_Subset(data_subset)` is assigned to `data`, but `FindClusters()` and the saved object are then applied to `data_subset`. This was preserved rather than silently corrected because changing it could alter the reproduced analysis. The eGFP threshold remains exactly `eGFP > 1.5`.

2. **`3.0_Subset_eGFP.R` also mixes full-object and subset visualizations.** Several downstream `FeaturePlot()`/`VlnPlot()` calls use `data` rather than `data_subset` or `data.harmony`. These calls were retained because they are part of the historical script.

3. **`3.2_Annotations.R` does not save the annotated Seurat object.** The original `saveRDS(data, ... "Data_Anotado.rds")` line is commented out, so this script only generates annotation figures unless that line is deliberately re-enabled.

4. **A Seurat-to-Python export step is not included among the supplied scripts.** `3.3_Generate_AnnData.py` expects `counts_<sample>.mtx`, `metadata_<sample>.csv/.xlsx`, `gene_names_<sample>.csv`, and `pca_<sample>.csv` under `Conversion_Seurat_To_Python/`. The script that creates those files was not supplied here.

5. **`3.4_Add_Velocity.py` depends on historical scVelo source objects.** Spliced, unspliced and ambiguous layers are transferred from `0.3_Dowstream_Analysis/TrayectoryAnalysis/AnnDataObj/AnnData_<sample>_scVelo.h5ad`. Those source objects are external dependencies of this folder.

6. **Barcode harmonization in `3.4_Add_Velocity.py` is potentially consequential.** The historical code truncates each cleaned barcode at the first hyphen (`bc.split('-')[0]`) before matching and then makes names unique. This behavior is preserved and should be verified against the original barcode conventions.

7. **`3.5_Velocity_New.py` performs two layers of preprocessing.** It first filters/normalizes/log-transforms/scales the combined AnnData object and constructs a PCA/neighborhood graph, then runs `scv.pp.filter_and_normalize()` separately in WT and MUT before dynamic velocity estimation. This historical analytical sequence was not altered.

8. **CellRank use in `3.5_Velocity_New.py` is limited to a `VelocityKernel`.** The script computes and stores a velocity-based transition matrix; it does not infer terminal states, fate probabilities or lineage drivers in the supplied version.

9. **`3.6_scVelo_Plots.py` contains multiple successive plotting attempts.** Earlier PDF/PNG variants (`Modify`, `V2`, `V3`) and later high-resolution `V10` variants are all retained for traceability. Before the final public release, decide which outputs correspond to the figures actually used in the manuscript/thesis and optionally remove obsolete plotting attempts.

10. **Utilities are still required.** The R workflow depends on `utils/0.3_Util_SeuratPipeline.R`, including functions such as `SeuratPipeline_Subset()` and `Calculo.Dimensiones.PCA()`. That utility was not part of the files supplied in this batch and therefore was not generated or reconstructed.

11. **Paths are now portable.** All scripts use the environment variable `PROJECT_ROOT` when available and otherwise default to the current working directory. The historical analytical subdirectory names were retained so dependencies remain traceable.
