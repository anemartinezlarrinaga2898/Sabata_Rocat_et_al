# Sabata_Rocat_et_al

**Context-dependent response of endothelial cells to PIK3CA mutation**

Welcome to this repository, where we delve into the code used to reproduce the analysis and generate all necessary R objects for the study: *"Endothelial Cell Fate and Differentiation Dictate Pathogenic Responses to PIK3CA Mutations in PROS"*.

**Abstract**

Cancer mutations in the PIK3CA gene cause congenital disorders. The endothelium is frequently affected in these disorders, displaying vascular overgrowth. Pathological PIK3CA vascular phenotypes are found in veins and capillaries but rarely in arteries for reasons that are unclear. Here, using lineage tracing, we show that expression of mutated PIK3CAH1047R in endothelial cells leads to clonal expansions in capillary and venous endothelial cells. In contrast, arterial endothelial cells are refractory to PIK3CAH1047R and never display pathological phenotypes. Moreover, PIK3CAH1047R expression in arterial precursors interrupts differentiation and drives fate switch towards venous identity. The PIK3CAH1047R-driven arterial-to-venous fate switch is orchestrated by upregulation of the vein-specifying transcription factor Nr2f2/COUP-TFII. Our findings reveal that pathogenic responses to PIK3CAH1047R greatly depend on the diferentation stage and fate trajectory of the targeted cell. Arteries are thus shielded against PIK3CA mutation, solving the long-standing question on the rarity of PIK3CA-related arterial malformations observed in patients.

## In this repository you would find: 

- **0.0_SeuratObjectGeneration.R** --> creation of the seurat object for each of the samples
- **0.1_Doublets.R** --> Doublet identificacion using DoubletFinder
- **0.2_AddMetadata_Doublets.R** --> Add the doublet information to the metadata of each samples
- **0.3_SeuratPipeline.R** --> Seurat Pipeline from normalization to clustering analysis + marker estimation
- **0.4_Clustering** --> Check for different resolutions in order to verify which fits better the data
- **0.5_RemoveClusters.R** --> Delete from the object the PECAM1- cells
- **0.6_Annotations.R**--> Generate the final annotations of the data
- **Plots_Paper.rmd** --> the different plots available in paper

### Computational approach summary: 

1. Aligment and quantification using mm10 mouse refenrence with CellRanger v 7.0.1
2. Pre-processing of the data to remove
  - Low quality cells and genes with no reads: nFeature >250, percent.mt >20, complexity >0.8
  - Doublets estimation using Doublet Finder (v.2.0.3)
3. Normalization, choosing highly variable gene, dimmensionality reduction
4. Clustering and annotation

**Data available at GEO: XXXXX**
