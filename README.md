# Sabata_Rocat_et_al

**Differential response of endothelial cells to PIK3CA mutations dictates vascular pathological phenotypes**

**Abstract**

Activating somatic mutations in PIK3CA are frequently found in cancer and in a group of congenital disorders known as PIK3CA-Related Overgrowth Spectrum (PROS). PIK3CA mutations induce large clone expansions that break tissue architecture. A predominant pattern of pathogenicity affecting mesodermal-derived tissues is observed but it is not understood. The vasculature is a tissue commonly affected in PROS that also shows context-dependent manifestation in certain endothelial cell types. To study why PIK3CA pathogenic phenotypes are restricted to some endothelial cell types, we combine a suite of inducible Cre lines with distinct Cre-dependent models expressing PIK3CAH1047R. We show that PIK3CAH1047R mutant endothelial cells overexpand in capillaries and veins but not in arteries. Expression of PIK3CAH1047R in the pre-arterial population, instead, ignites a switch towards a venous phenotype. This transition is orchestrated through the upregulation of vein-specifying transcription factor Nr2f2/COUP-TFII, which prevents arterial differentiation resulting in vascular malformations in capillaries and veins. Our findings underscore the importance of cell differentiation and cell fate in dictating the pathogenic response to PIK3CAH1047R mutations and solve a long-standing question regarding the rarity of PIK3CA-related arterial malformations in PROS patients. 

# In this repository you would find: 

- 0.0_SeuratObjectGeneration.R --> creation of the seurat object for each of the samples
- 0.1_Doublets.R --> Doublet identificacion using DoubletFinder
- 0.2_AddMetadata_Doublets.R --> Add the doublet information to the metadata of each samples
- 0.3_SeuratPipeline.R --> Seurat Pipeline from normalization to clustering analysis + marker estimation
- 0.4_Clustering --> Check for different resolutions in order to verify which fits better the data
- 0.5_RemoveClusters.R --> Delete from the object the PECAM1- cells
- 0.6_Annotations.R--> Generate the final annotations of the data
- Plots_Paper.rmd --> the different plots available in paper
