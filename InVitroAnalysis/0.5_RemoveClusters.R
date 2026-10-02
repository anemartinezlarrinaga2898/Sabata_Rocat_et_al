################################################################################
# SCRIPT: Removal of non-endothelial and low-quality clusters and reclustering
# AUTHOR: Ane Martinez Larrinaga
# DATE: 10-07-2024
#
# DESCRIPTION:
# This script removes clusters identified as non-endothelial cells and/or
# low-quality cell populations from the initial scRNA-seq clustering.
# Clusters 6 and 7 from the RNA_snn_res.0.3 clustering were excluded.
# Following their removal, variable feature selection, scaling, PCA,
# neighborhood graph construction, UMAP embedding and clustering were
# recalculated using the retained cells.
#
# Expression of established endothelial and vascular subtype markers was
# visualized to evaluate the resulting clusters. Cluster-specific markers
# were subsequently identified at resolutions 0.1 and 0.3.
#
# INPUT:
#   0.2_SeuratPipeline/Seu.Obj.rds
#
# OUTPUT:
#   0.2_SeuratPipeline/Seu.Obj_Remove.rds
#   UMAP plots after cluster removal
#   Sample and phenotype distribution plots
#   Endothelial marker FeaturePlots and VlnPlots
#   Cluster marker tables for resolutions 0.1 and 0.3
#
# CLUSTERS REMOVED:
#   RNA_snn_res.0.3 clusters 6 and 7
#   Reason: non-endothelial identity and/or low-quality transcriptional profile
#
# MAIN PARAMETERS:
#   Variable features: Seurat VST method
#   PCA dimensions: data-driven selection
#   Clustering resolutions: 0.1, 0.3 and 0.5
#   Marker detection: positive markers only
#   Minimum fraction of expressing cells: 0.25
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(Matrix)
library(matchSCore2)
library(writexl)
library(openxlsx)


# ------------------------------------------------------------------------------
# 2. Define output directory
# ------------------------------------------------------------------------------

path.guardar <- "0.2_SeuratPipeline"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}

getPalette <- colorRampPalette(
  brewer.pal(8, "Set1")
)


# ------------------------------------------------------------------------------
# 3. Load initial clustered Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  file.path(path.guardar, "Seu.Obj.rds")
)

data <- SetIdent(
  data,
  value = "RNA_snn_res.0.3"
)


# ------------------------------------------------------------------------------
# 4. Remove non-endothelial / low-quality clusters
# ------------------------------------------------------------------------------

# Clusters 6 and 7 were excluded based on their non-endothelial identity
# and/or low-quality transcriptional profiles.

data_cluster <- subset(
  data,
  idents = c("6", "7"),
  invert = TRUE
)

print(
  paste(
    "Cells retained after cluster removal:",
    ncol(data_cluster)
  )
)


# ------------------------------------------------------------------------------
# 5. Recalculate variable features
# ------------------------------------------------------------------------------

gene.info.distribution <- summary(
  Matrix::colSums(
    data_cluster@assays$RNA@counts > 0
  )
)

hvg.number <- round(
  gene.info.distribution[4] + 100
)

data_cluster <- FindVariableFeatures(
  object = data_cluster,
  selection.method = "vst",
  nfeatures = hvg.number
)

data_cluster <- ScaleData(
  object = data_cluster
)


# ------------------------------------------------------------------------------
# 6. Recalculate PCA
# ------------------------------------------------------------------------------

data_cluster <- RunPCA(
  object = data_cluster
)


# ------------------------------------------------------------------------------
# 7. Determine number of principal components
# ------------------------------------------------------------------------------

pct <- data_cluster[["pca"]]@stdev /
  sum(data_cluster[["pca"]]@stdev) * 100

cumu <- cumsum(pct)

co1 <- which(
  cumu > 90 & pct < 5
)[1]

co2 <- sort(
  which(
    (
      pct[1:(length(pct) - 1)] -
      pct[2:length(pct)]
    ) > 0.1
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
# 8. Recalculate neighbors, UMAP and clustering
# ------------------------------------------------------------------------------

data_cluster <- FindNeighbors(
  object = data_cluster,
  dims = 1:dim.final
)

data_cluster <- RunUMAP(
  object = data_cluster,
  dims = 1:dim.final
)

resolutions <- c(
  0.1,
  0.3,
  0.5
)

data_cluster <- FindClusters(
  object = data_cluster,
  resolution = resolutions
)


# ------------------------------------------------------------------------------
# 9. Visualize clustering after cluster removal
# ------------------------------------------------------------------------------

cluster.cols <- getPalette(
  length(unique(data_cluster$RNA_snn_res.0.5))
)

for (res in c("0.1", "0.3", "0.5")) {

  ident <- paste0(
    "RNA_snn_res.",
    res
  )

  p <- DimPlot(
    data_cluster,
    reduction = "umap",
    group.by = ident,
    label = TRUE,
    label.size = 5,
    cols = cluster.cols,
    pt.size = 1,
    raster = FALSE
  ) &
    NoAxes()

  if (res != "0.1") {
    p <- p & NoLegend()
  }

  ggsave(
    filename = file.path(
      path.guardar,
      paste0(ident, "_Remove.png")
    ),
    plot = p,
    width = 10,
    height = 10
  )
}


# ------------------------------------------------------------------------------
# 10. Visualize sample distribution
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
  data_cluster,
  reduction = "umap",
  group.by = "ID",
  label = FALSE,
  cols = scales::alpha(sample.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "ID_Remove.png"
  ),
  plot = p,
  width = 10,
  height = 10
)


# ------------------------------------------------------------------------------
# 11. Visualize phenotype distribution
# ------------------------------------------------------------------------------

phenotype.cols <- c(
  "#99C5E3",
  "#8CCE7D",
  "#FFC685"
)

p <- DimPlot(
  data_cluster,
  reduction = "umap",
  group.by = "Phenotype",
  split.by = "Phenotype",
  label = FALSE,
  cols = scales::alpha(phenotype.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "Phenotype_split_Remove.png"
  ),
  plot = p,
  width = 15,
  height = 7
)


# ------------------------------------------------------------------------------
# 12. Cluster composition by phenotype
# ------------------------------------------------------------------------------

composition.cols <- getPalette(10)

p <- matchSCore2::summary_barplot(
  class.fac = data_cluster$RNA_snn_res.0.1,
  obs.fac = data_cluster$Phenotype
) +
  scale_fill_manual(values = composition.cols)

ggsave(
  filename = file.path(
    path.guardar,
    "BarPlot_Pheno_Res01_Remove.png"
  ),
  plot = p,
  width = 5,
  height = 7
)


p <- matchSCore2::summary_barplot(
  class.fac = data_cluster$RNA_snn_res.0.3,
  obs.fac = data_cluster$Phenotype
) +
  scale_fill_manual(values = composition.cols)

ggsave(
  filename = file.path(
    path.guardar,
    "BarPlot_Pheno_Res03_Remove.png"
  ),
  plot = p,
  width = 5,
  height = 7
)


# ------------------------------------------------------------------------------
# 13. Visualize endothelial and vascular subtype markers
# ------------------------------------------------------------------------------

marker.groups <- list(

  Endothelial = c(
    "Pecam1", "Cdh5", "Vwf"
  ),

  Capillary = c(
    "Kdr", "Rgcc", "Cd200",
    "Cd300lg", "Cd36", "Sgk1"
  ),

  Arterial = c(
    "Sox17", "Hey1", "Sema3g", "Clu"
  ),

  Venous = c(
    "Nr2f2", "Vcam1", "Vwf", "Icam1"
  ),

  Angiogenic = c(
    "Esm1", "Cxcr4", "Dll4",
    "Col4a1", "Col4a2"
  ),

  Proliferative = c(
    "Mki67", "Cdk1", "Cdk2", "Cdk4", "Cdk6"
  ),

  Lymphatic = c(
    "Lyve1", "Prox1", "Pdpln"
  )
)


for (marker.name in names(marker.groups)) {

  genes <- marker.groups[[marker.name]]

  p <- FeaturePlot(
    data_cluster,
    features = genes,
    min.cutoff = "q9",
    order = TRUE,
    raster = FALSE,
    cols = c("Grey", "Red")
  ) &
    NoAxes()

  ggsave(
    filename = file.path(
      path.guardar,
      paste0(
        "FeaturePlot_",
        marker.name,
        "_Remove.png"
      )
    ),
    plot = p,
    width = 10,
    height = 10
  )
}


# ------------------------------------------------------------------------------
# 14. Generate violin plots of vascular subtype markers
# ------------------------------------------------------------------------------

data_cluster <- SetIdent(
  data_cluster,
  value = "RNA_snn_res.0.3"
)

vln.groups <- marker.groups[
  c(
    "Capillary",
    "Arterial",
    "Venous",
    "Angiogenic",
    "Proliferative",
    "Lymphatic"
  )
]

vln.cols <- getPalette(10)


for (marker.name in names(vln.groups)) {

  genes <- vln.groups[[marker.name]]

  p <- VlnPlot(
    object = data_cluster,
    features = genes,
    cols = vln.cols,
    pt.size = 0.1,
    sort = TRUE
  )

  ggsave(
    filename = file.path(
      path.guardar,
      paste0(
        "Vln_",
        marker.name,
        "_Remove.png"
      )
    ),
    plot = p,
    width = 10,
    height = 10
  )
}


# ------------------------------------------------------------------------------
# 15. Save reclustered endothelial Seurat object
# ------------------------------------------------------------------------------

data_cluster <- SetIdent(
  data_cluster,
  value = "RNA_snn_res.0.1"
)

saveRDS(
  data_cluster,
  file = file.path(
    path.guardar,
    "Seu.Obj_Remove.rds"
  )
)


# ------------------------------------------------------------------------------
# 16. Identify markers - resolution 0.1
# ------------------------------------------------------------------------------

data_cluster <- SetIdent(
  data_cluster,
  value = "RNA_snn_res.0.1"
)

markers.res01 <- FindAllMarkers(
  data_cluster,
  only.pos = TRUE,
  min.pct = 0.25
)

writexl::write_xlsx(
  markers.res01,
  file.path(
    path.guardar,
    "TotalMarkers_Res01_Remove.xlsx"
  )
)

markers.cluster.res01 <- split(
  markers.res01,
  markers.res01$cluster
)

openxlsx::write.xlsx(
  markers.cluster.res01,
  file.path(
    path.guardar,
    "Markers_ByCluster_01_Remove.xlsx"
  )
)


# ------------------------------------------------------------------------------
# 17. Identify markers - resolution 0.3
# ------------------------------------------------------------------------------

data_cluster <- SetIdent(
  data_cluster,
  value = "RNA_snn_res.0.3"
)

markers.res03 <- FindAllMarkers(
  data_cluster,
  only.pos = TRUE,
  min.pct = 0.25
)

writexl::write_xlsx(
  markers.res03,
  file.path(
    path.guardar,
    "TotalMarkers_Res03_Remove.xlsx"
  )
)

markers.cluster.res03 <- split(
  markers.res03,
  markers.res03$cluster
)

openxlsx::write.xlsx(
  markers.cluster.res03,
  file.path(
    path.guardar,
    "Markers_ByCluster_03_Remove.xlsx"
  )
)

print(
  "Removal of non-endothelial/low-quality clusters and reclustering completed."
)
