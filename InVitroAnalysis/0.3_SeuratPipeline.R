################################################################################
# SCRIPT: Seurat preprocessing, dimensionality reduction and clustering
# AUTHOR: Ane Martinez Larrinaga
# DATE: 23-10-2023
#
# DESCRIPTION:
# This script performs gene filtering, normalization, variable feature
# selection, scaling, PCA, dimensionality reduction and graph-based clustering
# of the scRNA-seq dataset. UMAP representations are generated to assess
# clustering, sample distribution, phenotype distribution and expression of
# selected endothelial cell markers. Cluster-specific marker genes are
# subsequently identified using Seurat.
#
# INPUT:
#   0.2_SeuratPipeline/Data_Doublets.rds
#
# OUTPUT:
#   0.2_SeuratPipeline/Seu.Obj.rds
#   UMAP and marker-expression plots
#   Cluster marker tables for resolutions 0.1 and 0.3
#
# MAIN PARAMETERS:
#   Genes retained if detected in > 5 cells
#   Variable features: Seurat VST method
#   PCA dimensions: data-driven selection
#   Clustering resolutions: 0.1, 0.3 and 0.5
#   Marker detection: positive markers, min.pct = 0.25
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(Matrix)
library(openxlsx)
library(writexl)
library(matchSCore2)

options(Seurat.object.assay.version = "v4")


# ------------------------------------------------------------------------------
# 2. Define output directory
# ------------------------------------------------------------------------------

path.guardar <- "0.2_SeuratPipeline"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 3. Load Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  file.path(path.guardar, "Data_Doublets.rds")
)


# ------------------------------------------------------------------------------
# 4. Filter lowly detected genes
# ------------------------------------------------------------------------------

# Retain genes detected in more than five cells
n_cells <- Matrix::rowSums(
  data@assays$RNA@counts > 0
)

kept_genes <- rownames(data)[
  n_cells > 5
]

data <- subset(
  data,
  features = kept_genes
)


# ------------------------------------------------------------------------------
# 5. Normalize data and identify variable features
# ------------------------------------------------------------------------------

data <- NormalizeData(
  object = data
)

# Determine the number of highly variable genes from the distribution
# of detected genes per cell
gene.info.distribution <- summary(
  Matrix::colSums(data@assays$RNA@counts > 0)
)

hvg.number <- round(
  gene.info.distribution[4] + 100
)

data <- FindVariableFeatures(
  object = data,
  selection.method = "vst",
  nfeatures = hvg.number
)

data <- ScaleData(
  object = data
)


# ------------------------------------------------------------------------------
# 6. Principal component analysis
# ------------------------------------------------------------------------------

data <- RunPCA(
  object = data
)


# ------------------------------------------------------------------------------
# 7. Determine number of principal components
# ------------------------------------------------------------------------------

pct <- data[["pca"]]@stdev /
  sum(data[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

# First PC at which cumulative explained variation exceeds 90%
# while individual contribution is below 5%
co1 <- which(
  cumu > 90 & pct < 5
)[1]

# Identify the last major drop in variation between consecutive PCs
co2 <- sort(
  which(
    (pct[1:(length(pct) - 1)] -
       pct[2:length(pct)]) > 0.1
  ),
  decreasing = TRUE
)[1] + 1

dim.final <- min(
  co1,
  co2
)

print(
  paste("Number of PCs selected:", dim.final)
)


# ------------------------------------------------------------------------------
# 8. Graph construction, UMAP and clustering
# ------------------------------------------------------------------------------

data <- FindNeighbors(
  object = data,
  dims = 1:dim.final
)

data <- RunUMAP(
  object = data,
  dims = 1:dim.final
)

resolutions <- c(
  0.1,
  0.3,
  0.5
)

data <- FindClusters(
  object = data,
  resolution = resolutions
)


# ------------------------------------------------------------------------------
# 9. Visualize clustering resolutions
# ------------------------------------------------------------------------------

cluster.palette <- colorRampPalette(
  brewer.pal(8, "Set1")
)

col <- cluster.palette(
  length(unique(data$RNA_snn_res.0.5))
)


p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "RNA_snn_res.0.1",
  label = TRUE,
  label.size = 5,
  cols = col,
  pt.size = 1,
  raster = FALSE
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "RNA_snn_res.0.1.png"),
  plot = p,
  width = 10,
  height = 10
)


p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "RNA_snn_res.0.3",
  label = TRUE,
  label.size = 5,
  cols = col,
  pt.size = 1,
  raster = FALSE
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(path.guardar, "RNA_snn_res.0.3.png"),
  plot = p,
  width = 10,
  height = 10
)


p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "RNA_snn_res.0.5",
  label = TRUE,
  label.size = 5,
  cols = col,
  pt.size = 1,
  raster = FALSE
) & NoAxes() & NoLegend()

ggsave(
  filename = file.path(path.guardar, "RNA_snn_res.0.5.png"),
  plot = p,
  width = 10,
  height = 10
)


# ------------------------------------------------------------------------------
# 10. Visualize sample and phenotype distribution
# ------------------------------------------------------------------------------

sample.cols <- c(
  "#99C5E3",
  "#8CCE7D",
  "#FFC685",
  "#E2B5D5",
  "#AAAEB0",
  "#FA8D76",
  "#F4D166"
)

p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "ID",
  label = FALSE,
  cols = scales::alpha(sample.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "ID.png"),
  plot = p,
  width = 10,
  height = 10
)


phenotype.cols <- c(
  "#99C5E3",
  "#8CCE7D",
  "#FFC685"
)

p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "Phenotype",
  label = FALSE,
  cols = scales::alpha(phenotype.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Phenotype.png"),
  plot = p,
  width = 10,
  height = 10
)


p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "Phenotype",
  split.by = "Phenotype",
  label = FALSE,
  cols = scales::alpha(phenotype.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Phenotype_split.png"),
  plot = p,
  width = 15,
  height = 7
)


# ------------------------------------------------------------------------------
# 11. Cluster composition by phenotype
# ------------------------------------------------------------------------------

composition.cols <- cluster.palette(10)

p <- matchSCore2::summary_barplot(
  class.fac = data$RNA_snn_res.0.1,
  obs.fac = data$Phenotype
) +
  scale_fill_manual(values = composition.cols)

ggsave(
  filename = file.path(path.guardar, "BarPlot_Pheno_Res01.png"),
  plot = p,
  width = 5,
  height = 7
)


p <- matchSCore2::summary_barplot(
  class.fac = data$RNA_snn_res.0.3,
  obs.fac = data$Phenotype
) +
  scale_fill_manual(values = composition.cols)

ggsave(
  filename = file.path(path.guardar, "BarPlot_Pheno_Res03.png"),
  plot = p,
  width = 5,
  height = 7
)


# ------------------------------------------------------------------------------
# 12. Expression of endothelial and subtype marker genes
# ------------------------------------------------------------------------------

# General endothelial markers
p <- FeaturePlot(
  data,
  features = c("Pecam1", "Cdh5", "Vwf"),
  min.cutoff = "q9",
  order = TRUE,
  pt.size = 1,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "FeaturePlots_MarkersEndo.png"),
  plot = p,
  width = 10,
  height = 10
)


# Capillary markers
p <- FeaturePlot(
  data,
  features = c(
    "Kdr", "Rgcc", "Cd200",
    "Cd300lg", "Cd36", "Sgk1"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Capillary.png"),
  plot = p,
  width = 10,
  height = 10
)


# Arterial markers
p <- FeaturePlot(
  data,
  features = c(
    "Sox17", "Hey1", "Sema3g", "Clu"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Artery.png"),
  plot = p,
  width = 10,
  height = 10
)


# Venous markers
p <- FeaturePlot(
  data,
  features = c(
    "Nr2f2", "Vcam1", "Vwf", "Icam1"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Venous.png"),
  plot = p,
  width = 10,
  height = 10
)


# Angiogenic / tip-cell markers
p <- FeaturePlot(
  data,
  features = c(
    "Esm1", "Cxcr4", "Dll4",
    "Col4a1", "Col4a2"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Tip.png"),
  plot = p,
  width = 10,
  height = 10
)


# Proliferation markers
p <- FeaturePlot(
  data,
  features = c(
    "Mki67", "Cdk1", "Cdk2", "Cdk4", "Cdk6"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Division.png"),
  plot = p,
  width = 10,
  height = 10
)


# Lymphatic endothelial markers
p <- FeaturePlot(
  data,
  features = c("Lyve1", "Prox1", "Pdpln"),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Lymphatics.png"),
  plot = p,
  width = 10,
  height = 5
)


# Additional endothelial markers of interest
p <- FeaturePlot(
  data,
  features = c(
    "Lyve1", "Prox1", "Hey1", "Hey2",
    "Unc5b", "Flt4", "Nrp2", "Nr2f2",
    "Ephb4", "Cxcr4"
  ),
  min.cutoff = "q9",
  order = TRUE,
  raster = FALSE,
  cols = c("Grey", "Red")
) & NoAxes()

ggsave(
  filename = file.path(path.guardar, "Fp_Interest.png"),
  plot = p,
  width = 15,
  height = 10
)


# ------------------------------------------------------------------------------
# 13. Save processed Seurat object
# ------------------------------------------------------------------------------

data <- SetIdent(
  data,
  value = "RNA_snn_res.0.1"
)

saveRDS(
  data,
  file = file.path(path.guardar, "Seu.Obj.rds")
)


# ------------------------------------------------------------------------------
# 14. Identify cluster markers - resolution 0.1
# ------------------------------------------------------------------------------

data <- SetIdent(
  data,
  value = "RNA_snn_res.0.1"
)

markers.res01 <- FindAllMarkers(
  data,
  only.pos = TRUE,
  min.pct = 0.25
)

writexl::write_xlsx(
  markers.res01,
  file.path(path.guardar, "TotalMarkers_Res01.xlsx")
)

markers.cluster.res01 <- split(
  markers.res01,
  markers.res01$cluster
)

openxlsx::write.xlsx(
  markers.cluster.res01,
  file.path(path.guardar, "Markers_ByCluster.xlsx")
)


# ------------------------------------------------------------------------------
# 15. Identify cluster markers - resolution 0.3
# ------------------------------------------------------------------------------

data <- SetIdent(
  data,
  value = "RNA_snn_res.0.3"
)

markers.res03 <- FindAllMarkers(
  data,
  only.pos = TRUE,
  min.pct = 0.25
)

writexl::write_xlsx(
  markers.res03,
  file.path(path.guardar, "TotalMarkers_Res03.xlsx")
)

markers.cluster.res03 <- split(
  markers.res03,
  markers.res03$cluster
)

openxlsx::write.xlsx(
  markers.cluster.res03,
  file.path(path.guardar, "Markers_ByCluster_03.xlsx")
)

print("Seurat preprocessing and clustering completed.")
