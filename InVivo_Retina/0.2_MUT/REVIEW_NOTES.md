# REVIEW NOTES — MUT preprocessing scripts

These notes document points that should be checked before treating the scripts as the final archival version of the analysis. The scripts were homogenized for readability and GitHub publication without intentionally changing analytical thresholds, clustering resolutions, selected clusters, marker sets, or biological decisions.

## 1. Project-relative paths

All absolute `/ijc/...` paths and the private `Paths.R` dependency were removed. The scripts now use:

```r
project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
results_root <- file.path(project_root, "results")
utils_dir <- file.path(project_root, "utils")
```

The utilities referenced by the original analysis therefore need to be available in `utils/`, including at least:

- `0.2_Util_DoubletDetection.R`
- `0.3_Util_SeuratPipeline.R`
- `0.19_Util_CellType_Classifier.R` where required.

## 2. MUT sample selection in `0.0_SeuratObject.R`

The original script extracts samples `L1152` and `L1155` from a previously filtered global retinal object. This selection has been preserved and made explicit as:

```r
mut_samples <- c("L1152", "L1155")
```

The input is expected at `data/Data.Filtered.rds`. Confirm that this is the intended public input object and that these sample identifiers are the ones you want exposed in the repository.

## 3. DoubletFinder predictions are retained

`0.1_Doublets.R` estimates doublets independently by sample and adds the DoubletFinder results to the metadata before merging the objects. The script does not remove predicted doublets. This behavior has been preserved.

## 4. `SeuratPipeline()` and `SeuratPipeline_Subset()` are external project functions

Several preprocessing decisions are encapsulated in custom utility functions. Their exact normalization, variable-feature selection, PCA and filtering behavior cannot be inferred from these scripts alone. The corresponding utility scripts must be included in the final repository if full reproducibility is required.

## 5. `1.0_Subset_EC.R`: plotting object after Harmony

The original script performs Harmony on `data_interest` and creates `data.harmony`, but some subsequent marker plots use `data` rather than `data.harmony`. This has been left unchanged because changing the plotted object would alter the historical script behavior. Confirm whether those plots were intended to show the pre-subset/pre-integration object or the Harmony-integrated endothelial subset.

## 6. `1.2_eGFP_Dist.R`: original path typo corrected

The original output path starts with `.5_MUT_Analysis/...` instead of `0.5_MUT_Analysis/...`. This was treated as an obvious path typo and standardized to `0.5_MUT_Analysis/...`; no analytical calculation is affected.

## 7. `1.2_eGFP_Dist.R` and `1.5_eGFP_Dist.R` contain a WT-specific exploratory section

Although both scripts are located in the MUT workflow, the original code explicitly performs:

```r
data <- SetIdent(data, value = "Phenotype")
data_wt <- subset(data, ident = "WT")
```

and then calculates ccAFv2 cell-cycle states, a pre-arterial UCell signature, and correlations in `data_wt`.

This section has intentionally been retained rather than silently changed to MUT. Confirm whether:
- this WT analysis is intentional and belongs in these scripts, or
- it is copied from a WT workflow and should be removed/replaced in the final repository.

## 8. `PredictCellCycle()` dependency

The eGFP scripts call `PredictCellCycle()`, apparently supplied by `0.19_Util_CellType_Classifier.R`. The utility must therefore be included or the function source documented.

## 9. eGFP-positive threshold

Both eGFP scripts define positivity using:

```r
FetchData(data, vars = "eGFP")$eGFP > 1.5
```

The threshold `1.5` has been preserved exactly. The factor comparisons were standardized to explicit string comparisons (`"TRUE"` / `"FALSE"`) for clarity without changing the intended classification.

## 10. Marker scripts require a command-line identity

`0.3_Marker.R`, `0.5_Markers.R`, `1.1_Marker.R`, `1.2_eGFP_Dist.R`, `1.4_Markers.R`, and `1.5_eGFP_Dist.R` depend on a metadata-column name supplied as a command-line argument. Input validation has been added, but the required identity for the published analysis should be documented when the final repository structure is assembled.

## 11. Historical analytical choices preserved

The following selections were preserved exactly from the supplied scripts:

- Complete MUT Harmony integration across `ID`.
- Initial endothelial subset from `Harmony_Log_res.0.3`: clusters `0, 1, 2, 3, 4, 9`.
- Subsequent exclusion from `Harmony_Log_res.0.5`: clusters `3, 6, 7, 9`.
- Original clustering-resolution vectors.
- Original marker gene panels and eGFP threshold.
- Original `FindMarkers` / `FindAllMarkers` calls and their explicitly supplied parameters.

No attempt was made to retrospectively optimize these choices.

## 12. `0.2_SeuratPipeline.R`: `cols_pheno` is not defined locally

The phenotype UMAP uses `cols_pheno`, but this object is not defined inside the supplied script. It may be created by one of the sourced utility files or by the original interactive environment. Confirm its source before expecting the script to run in a clean R session.

## 13. Doublet metadata alignment

`0.1_Doublets.R` adds the DoubletFinder result data frame to each sample using `cbind()` with the Seurat metadata. This assumes that the rows of the DoubletFinder result are in exactly the same order as the Seurat metadata. The historical implementation has been preserved. For a stricter reproducibility version, row-name alignment should be explicitly checked before combining the tables.

## 14. Active identity in the eGFP scripts

In `1.2_eGFP_Dist.R` and `1.5_eGFP_Dist.R`, some `summary_barplot()` calls use `data@active.ident` without first setting the active identity to the command-line `ident` variable. The plots may therefore depend on the active identity stored in the input RDS rather than the identity passed to the script. This has not been silently changed.

## 15. Hard-coded identity-marker panel in the eGFP scripts

The eGFP scripts contain a manually defined marker panel described by cluster numbers 0–5. Because the script accepts `ident` as a command-line argument, confirm that the marker panel corresponds to the same clustering resolution used for the published figure.

## 16. Marker name to verify

The supplied vascular-marker panels contain the string `Pdpln`. This was preserved exactly. Confirm that this is the intended gene symbol before final publication.
