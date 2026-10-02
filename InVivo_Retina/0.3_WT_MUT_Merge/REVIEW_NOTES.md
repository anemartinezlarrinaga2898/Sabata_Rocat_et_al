# REVIEW NOTES — combined WT + MUT scripts

These notes document points that should be checked before treating the scripts as the final archival version of the combined in vivo WT + MUT analysis. The files were homogenized for readability and GitHub publication without intentionally changing biological selections, thresholds, cluster mappings or statistical parameters.

## 1. Project-relative paths

All private `/ijc/...` paths and the `Paths.R` dependency were removed. Scripts now assume the repository is run from its root (or that `PROJECT_ROOT` is set) and use:

```r
project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
results_root <- file.path(project_root, "results")
utils_dir <- file.path(project_root, "utils")
data_dir <- file.path(project_root, "data")
```

The original utility `0.3_Util_SeuratPipeline.R` must therefore be available under `utils/` for scripts that call custom Seurat helper functions.

## 2. `0.1_Join_Obj.R`: WT and MUT inputs are not symmetric

The supplied script merges:

- a final annotated WT object (`Data_WT_Anotado.rds`), and
- the MUT Harmony object (`Harmony.rds`).

This asymmetry was preserved because it reflects the supplied historical workflow. Confirm that these are the exact objects used to generate the published joint analysis.

## 3. eGFP-positive threshold

Where eGFP positivity is calculated, the original threshold is preserved exactly:

```r
eGFP > 1.5
```

No attempt was made to re-estimate or optimize this threshold.

## 4. `0.1_Join_Obj.R`: inherited WT annotation metadata

After WT and MUT are merged and Harmony is rerun, the script visualizes `AnnotLayer`. Because this metadata originated from the WT object, MUT cells may initially have missing annotation values. The original handling of missing values as `No_Annot` was retained.

The script also saves several `AnnotLayer_WT.png` plots to the same filename, so later plots overwrite earlier versions. This historical behavior was not changed.

## 5. `0.3_eGFP_Levels.R`: Seurat syntax normalized

The supplied script used:

```r
SetIdent(data, ident = "eGFP_pos")
subset(data, value = "TRUE")
```

These argument names are not the standard Seurat interface. They were normalized to:

```r
SetIdent(data, value = "eGFP_pos")
subset(data, idents = "TRUE")
```

This is intended to implement the apparent original operation (select eGFP-positive cells) without changing the biological threshold.

## 6. Curated external inputs

Two absolute-path Excel inputs were replaced by repository-relative inputs under `data/`:

- `Annot_PCA_24_Info.xlsx` (sheet `markers`) used by `0.4_DotPlot_Annots.R`.
- `WT_MUT_markers.xlsx` (sheet `Hoja2`) used by `4.0_Annotate.R`.

These exact curated files need to be included in the final repository or documented as supplementary inputs.

## 7. `4.0_Annotate.R` and `4.2_NewAnnotation.R` are different annotation schemes

The two scripts apply different manual cluster mappings. `4.0_Annotate.R` writes to `EC_Subset/Annotations`, whereas `4.2_NewAnnotation.R` writes to `EC_Subset/New_Annotations`.

This distinction was intentionally preserved. Confirm which annotation version corresponds to each final figure/analysis and whether both should remain in the public repository.

## 8. Downstream scripts use different annotation branches

The supplied downstream scripts do not all consume the same annotation object:

- `5.0_DEG_Analysis.R` reads `EC_Subset/New_Annotations/Data_EC_Annotated.rds`.
- `4.3_Heatmap.R`, `5.4_Downsampling.R`, and `5.5_Pseudobulk.R` read `EC_Subset/Annotations/Data_EC_Annotated.rds`.
- `5.2_DotPlot_DEG.R` reads results from `DownStream/DEG/Annotations`, while `5.0_DEG_Analysis.R` writes to `DownStream/DEG/New_Annotations`.

No branch was silently redirected. This should be reconciled against the final figures and Methods before publication.

## 9. `4.3_Heatmap.R`: `colorRamp2` dependency corrected

The original script attempted to load `library(colorRamp2)`, although `colorRamp2()` is provided by the `circlize` package. The import was standardized to `library(circlize)`. The heatmap calculations and parameters were not changed.

## 10. `5.2_DotPlot_DEG.R`: `mat2` is calculated but not plotted

The script creates a winsorized matrix:

```r
mat2 <- pmax(pmin(mat, 2), -2)
```

but the heatmap is generated from `mat`, not `mat2`. This was left unchanged because switching the plotted matrix would alter the historical output. Confirm whether `mat2` was intended for the final figure.

`ComplexHeatmap` was added as an explicit library because the script calls `Heatmap()`.

## 11. `5.0_DEG_Analysis.R`: comparison being performed

Within each `Layer_1` population, differential expression is calculated only when both `MUT_TRUE` and `WT_TRUE` contain at least three cells. The contrast is:

```r
FindMarkers(
  data_c,
  ident.1 = "MUT_TRUE",
  ident.2 = "WT_TRUE"
)
```

Thus this script compares MUT versus WT specifically among cells classified as eGFP-positive. This logic was preserved.

## 12. `5.4_Downsampling.R`: 1,000 single-cell resampling rounds

The supplied analysis balances the number of WT and MUT eGFP-positive cells to the smaller group within a selected endothelial population, then repeats `FindMarkers()` 1,000 times with:

```r
logfc.threshold = 0
min.pct = 0
```

The random seed (`123`) and all analytical parameters were preserved. Only command-line argument validation was added.

## 13. `5.5_Pseudobulk.R`: singlet filtering differs from earlier single-cell analyses

The pseudobulk script explicitly subsets to cells classified as `Singlet` before aggregation. This differs from the earlier preprocessing pipelines, where DoubletFinder classifications were retained without removing predicted doublets. The pseudobulk-specific filtering was preserved exactly and should be described consistently if this analysis is reported.

## 14. `5.5_Pseudobulk.R`: sample-name parsing should be verified

The script constructs sample labels from `eGFP_pheno` and `ID`, aggregates counts by `cell` and `samples`, and subsequently parses aggregated column names using regular expressions before searching for `TRUE`, `MUT` and `WT` labels.

Because the exact column-name format produced by `AggregateExpression()` can depend on the grouping labels and Seurat behavior/version, verify that the generated `colData` contains the intended biological samples and conditions before treating the pseudobulk output as final.

## 15. `5.5_Pseudobulk.R`: unused parameters

`min_cells` and `min_pct` are defined at the beginning of the script but are not subsequently used in the supplied code. They were retained to avoid retroactively changing the analysis.

## 16. GSEA output path

The final line of `5.5_Pseudobulk.R` writes `GSEA_selected_top15.xlsx` without `path.guardar`, unlike the other pseudobulk outputs. This was preserved rather than silently relocating a historical output. Confirm whether it should instead be written inside the cell-type pseudobulk results directory.

## 17. Manual biological annotations were not altered

Cluster-to-cell-type mappings, endothelial labels, marker lists, eGFP threshold, palette choices, DE cutoffs, downsampling rounds, DESeq2 design, enrichment thresholds and GSEA filters were preserved as supplied. The scripts were cleaned for structure, portability and readability rather than re-analysed or optimized.
