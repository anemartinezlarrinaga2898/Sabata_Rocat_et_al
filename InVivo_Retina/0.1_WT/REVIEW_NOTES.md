# Review notes before public release

## Changes made in this GitHub-ready copy

1. Removed institution/server-specific `/ijc/...` paths and the private `Paths.R` dependency; paths now come from `R/project_config.R`.
2. `0.0_SeuratObject.R`: normalized its output to the WT object directory consumed by `0.1_Doublets.R`. The supplied file pointed to `0.5_MUT_Analysis/Obj`; confirm that this uploaded `0.0` file and `SAMPLES_TO_KEEP = c("L1152", "L1155")` are the intended WT first step.
3. `1.2_Clean_EC_Comparment.R`: normalized the Harmony output path to `Res03_Clus67/Harmony`, which is the path consumed by `1.3_Markers.R`. No analytical parameter changed.
4. `1.0_Subset_EC.R`: the final endothelial marker plots now use the integrated EC object (`data.harmony`) rather than the pre-subset whole object. This affects plotting only.
5. `3.3_Remove_Sample.R`: removed a duplicated `1.5` entry from the clustering resolution vector; the set of tested resolutions is unchanged.
6. `5.0_Annotate_WT.R`: fixed palette-name assignment and a `ggsave(height=...)` typo.
7. `5.1_Heatmap.R`: replaced `library(colorRamp2)` with `library(circlize)` and explicit `circlize::colorRamp2()`.
8. `5.2_Similarity.R`: corrected a plot title that referred to pericytes/tissues; calculations were not changed.

## Missing files / intermediate steps

The exact original versions of these utility scripts are required and were not supplied:
- `R/0.2_Util_DoubletDetection.R`
- `R/0.3_Util_SeuratPipeline.R`
- `R/0.19_Util_CellType_Classifier.R`
- `R/0.5_Util_Annotations.R`

The supplied numbered files also reference intermediate objects whose generating scripts were not included:
- `Data_Prolif_Signature_Estimated.rds` before `3.0_remove_New_Clusters.R`
- `Test_PCA_Dimns/PCA_Dimns_24/Harmony.rds` before `5.0_Annotate_WT.R`

Curated inputs also need to be added:
- `config/WT_Annot.xlsx`
- `config/WT_Final_Annot.xlsx`
- `config/sample_metadata.csv`

## Confirm before release

- Confirm sample IDs/phenotypes in `0.0_SeuratObject.R`.
- Confirm all manual cluster/sample exclusions match the final Methods and figures.
- Confirm the `eGFP > 1.5` threshold is documented consistently in the manuscript.
- Add exact R/package versions from the actual analysis environment.
