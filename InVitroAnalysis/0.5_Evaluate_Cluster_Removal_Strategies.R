################################################################################
# SCRIPT: Evaluation of cluster removal strategies
# AUTHOR: Ane Martinez Larrinaga
# DATE: 10-07-2024
#
# DESCRIPTION:
# This script evaluates alternative strategies for removing clusters considered
# non-endothelial, low-quality or transcriptionally ambiguous from the initial
# scRNA-seq dataset.
#
# Starting from clustering resolution 0.8, three alternative cluster-removal
# strategies are evaluated:
#
#   Option 1: clusters 7 and 9
#   Option 2: clusters 7, 9 and 8
#   Option 3: clusters 5, 7, 9 and 8
#
# Following each exclusion strategy, variable feature selection, scaling, PCA,
# neighborhood graph construction, UMAP embedding and clustering are
# recalculated.
#
# The resulting cell populations are evaluated using UCell scores derived from
# endothelial and vascular subtype marker signatures. Cell-cycle scores are
# additionally calculated to assess the contribution of proliferative states.
#
# INPUT:
#   0.2_SeuratPipeline/Seu.Obj.rds
#
# OUTPUT:
#   Results for each cluster-removal strategy, including:
#     - UCell signature FeaturePlots
#     - UCell signature violin plots
#     - Cell-cycle phase UMAP
#     - S-phase score FeaturePlot
#     - G2/M-phase score FeaturePlot
#
# CLUSTERING RESOLUTIONS:
#   0.1, 0.3, 0.5, 0.7, 0.8 and 0.9
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(Matrix)
library(gridExtra)
library(patchwork)
library(UCell)


# ------------------------------------------------------------------------------
# 2. Define paths and plotting parameters
# ------------------------------------------------------------------------------

base.output <- "0.2_SeuratPipeline/Clustering/RemovingCluster"

if (!dir.exists(base.output)) {
  dir.create(base.output, recursive = TRUE)
}

getPalette <- colorRampPalette(
  brewer.pal(8, "Set1")
)

col <- getPalette(20)


# ------------------------------------------------------------------------------
# 3. Function for reclustering after cluster removal
# ------------------------------------------------------------------------------

SeuratPipeline_Subset <- function(data_cluster) {

  # Determine number of highly variable genes
  gene.info.distribution <- summary(
    Matrix::colSums(
      data_cluster@assays$RNA@counts > 0
    )
  )

  hvg.number <- round(
    gene.info.distribution[4] + 100
  )

  # Variable feature selection
  data_cluster <- FindVariableFeatures(
    object = data_cluster,
    selection.method = "vst",
    nfeatures = hvg.number
  )

  # Scaling
  data_cluster <- ScaleData(
    object = data_cluster
  )

  # PCA
  data_cluster <- RunPCA(
    object = data_cluster
  )

  # Estimate number of PCs
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

  # Neighbors and UMAP
  data_cluster <- FindNeighbors(
    object = data_cluster,
    dims = 1:dim.final
  )

  data_cluster <- RunUMAP(
    object = data_cluster,
    dims = 1:dim.final
  )

  # Clustering
  resolutions <- c(
    0.1,
    0.3,
    0.5,
    0.7,
    0.8,
    0.9
  )

  data_cluster <- FindClusters(
    object = data_cluster,
    resolution = resolutions
  )

  return(data_cluster)
}


# ------------------------------------------------------------------------------
# 4. Function for UCell-based annotation assessment
# ------------------------------------------------------------------------------

Annotation_Function <- function(
    data,
    marker.list,
    output.path) {

  features_df <- do.call(
    rbind,
    marker.list
  )

  cell.types <- names(
    marker.list
  )

  signature.list <- list()

  for (i in seq_along(marker.list)) {

    cell.type <- cell.types[i]

    genes <- features_df$Features[
      features_df$CellType == cell.type
    ]

    signature.list[[cell.type]] <- genes
  }

  # Calculate UCell signature scores
  data <- UCell::AddModuleScore_UCell(
    data,
    features = signature.list
  )

  signature.names <- paste0(
    names(signature.list),
    "_UCell"
  )

  feature.plots <- list()
  violin.plots <- list()

  for (i in seq_along(signature.names)) {

    signature <- signature.names[i]

    feature.plots[[i]] <- FeaturePlot(
      data,
      features = signature,
      min.cutoff = "q9",
      order = TRUE
    ) &
      NoAxes()

    violin.plots[[i]] <- VlnPlot(
      data,
      features = signature,
      sort = TRUE,
      cols = col
    ) &
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      )
  }

  p <- do.call(
    "grid.arrange",
    c(
      feature.plots,
      ncol = 5,
      nrow = 2
    )
  )

  ggsave(
    filename = file.path(
      output.path,
      "FeaturePlot_Annotation.png"
    ),
    plot = p,
    width = 30,
    height = 12
  )

  p <- do.call(
    "grid.arrange",
    c(
      violin.plots,
      ncol = 5,
      nrow = 2
    )
  )

  ggsave(
    filename = file.path(
      output.path,
      "ViolinPlot_Annotation.png"
    ),
    plot = p,
    width = 30,
    height = 12
  )

  return(data)
}


# ------------------------------------------------------------------------------
# 5. Load Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  "0.2_SeuratPipeline/Seu.Obj.rds"
)


# ------------------------------------------------------------------------------
# 6. Recalculate clustering used for cluster-removal evaluation
# ------------------------------------------------------------------------------

# NOTE:
# Calculo.Dimensiones.PCA() must be available in the accompanying utility
# scripts included in the repository.

dim.final <- Calculo.Dimensiones.PCA(
  data
)

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
  0.5,
  0.7,
  0.8,
  0.9
)

data <- FindClusters(
  data,
  resolution = resolutions
)

# Resolution 0.8 was used to evaluate clusters for exclusion
data <- SetIdent(
  data,
  value = "RNA_snn_res.0.8"
)


# ------------------------------------------------------------------------------
# 7. Define endothelial and vascular marker signatures
# ------------------------------------------------------------------------------

EndothelialCells <- data.frame(
  Features = c("Pecam1", "Cdh5"),
  CellType = "EndothelialCells"
)

Artery <- data.frame(
  Features = c(
    "Sox17", "Dkk2", "Efbn2", "Sema3g",
    "Fbln5", "Hey1", "Bmx", "Mgp"
  ),
  CellType = "Artery"
)

Venous <- data.frame(
  Features = c(
    "Ackr1", "Selp", "Sele", "Icam1",
    "Vwf", "Vcam1", "Nr2f2", "Ephb4"
  ),
  CellType = "Venous"
)

Capillary <- data.frame(
  Features = c(
    "Rgcc", "Cd36", "Kdr", "Car2",
    "Cd200", "Cd300lg", "Giphbp1", "Aqp1"
  ),
  CellType = "Capillary"
)

TipCells <- data.frame(
  Features = c(
    "Esm1", "Dll4", "Cxcr4", "Col4a2",
    "Col4a1", "Sparc", "Insr", "Angpt2"
  ),
  CellType = "TipCells"
)

Lymphatics <- data.frame(
  Features = c(
    "Prox1", "Lyve1", "Pdpn", "Cd63",
    "Flt4", "Ccdc3", "Mmrn2", "Foxp2",
    "Cldn11", "Net1", "Neo1"
  ),
  CellType = "Lymphatics"
)

Fenestration <- data.frame(
  Features = c("Plvap", "Plpp3", "Igfbp3"),
  CellType = "Fenestration"
)

Angiogenesis <- data.frame(
  Features = c(
    "Angpt1", "Pdgfc", "Vegfc", "Vegfb",
    "Lrp1", "Tek", "Kdr", "Angpt2"
  ),
  CellType = "Angiogenesis"
)

Immature <- data.frame(
  Features = c(
    "Pdlim1", "Rplp0", "Gapdh", "Igfbp7",
    "Rps2", "Rpsa", "Rpl12", "Aplnr"
  ),
  CellType = "Immature"
)

Lista_Markers <- list(
  EndothelialCells = EndothelialCells,
  Artery = Artery,
  Venous = Venous,
  Capillary = Capillary,
  TipCells = TipCells,
  Lymphatics = Lymphatics,
  Fenestration = Fenestration,
  Angiogenesis = Angiogenesis,
  Immature = Immature
)


# ------------------------------------------------------------------------------
# 8. Prepare cell-cycle genes
# ------------------------------------------------------------------------------

s.genes <- stringr::str_to_title(
  cc.genes$s.genes
)

g2m.genes <- stringr::str_to_title(
  cc.genes$g2m.genes
)


# ------------------------------------------------------------------------------
# 9. Define alternative cluster-removal strategies
# ------------------------------------------------------------------------------

removal.strategies <- list(

  Remove_Clus_7_9 =
    c("7", "9"),

  Remove_Clus_7_9_8 =
    c("7", "9", "8"),

  Remove_Clus_7_9_8_5 =
    c("7", "9", "8", "5")
)


# ------------------------------------------------------------------------------
# 10. Evaluate each cluster-removal strategy
# ------------------------------------------------------------------------------

for (strategy.name in names(removal.strategies)) {

  clusters.to.remove <-
    removal.strategies[[strategy.name]]

  print(
    paste(
      "Evaluating strategy:",
      strategy.name
    )
  )

  print(
    paste(
      "Removing clusters:",
      paste(
        clusters.to.remove,
        collapse = ", "
      )
    )
  )


  # --------------------------------------------------------------------------
  # Remove selected clusters
  # --------------------------------------------------------------------------

  data_cluster <- subset(
    data,
    idents = clusters.to.remove,
    invert = TRUE
  )


  # --------------------------------------------------------------------------
  # Define strategy-specific output directory
  # --------------------------------------------------------------------------

  output.path <- file.path(
    base.output,
    strategy.name
  )

  if (!dir.exists(output.path)) {
    dir.create(
      output.path,
      recursive = TRUE
    )
  }


  # --------------------------------------------------------------------------
  # Reprocess and recluster retained cells
  # --------------------------------------------------------------------------

  data_cluster <- SeuratPipeline_Subset(
    data_cluster
  )


  # --------------------------------------------------------------------------
  # Evaluate endothelial signatures
  # --------------------------------------------------------------------------

  data_cluster <- Annotation_Function(
    data_cluster,
    Lista_Markers,
    output.path
  )


  # --------------------------------------------------------------------------
  # Cell-cycle scoring
  # --------------------------------------------------------------------------

  data_cluster <- CellCycleScoring(
    data_cluster,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )


  # --------------------------------------------------------------------------
  # Cell-cycle phase UMAP
  # --------------------------------------------------------------------------

  p <- DimPlot(
    data_cluster,
    reduction = "umap",
    group.by = "Phase",
    label = TRUE,
    label.size = 5,
    cols = c(
      "#95BDED",
      "#F5A040",
      "#C8C2F2"
    ),
    pt.size = 0.5,
    raster = FALSE
  ) &
    NoAxes()

  ggsave(
    filename = file.path(
      output.path,
      "Phase.png"
    ),
    plot = p,
    width = 7,
    height = 7
  )


  # --------------------------------------------------------------------------
  # S-phase score
  # --------------------------------------------------------------------------

  p <- FeaturePlot(
    data_cluster,
    features = "S.Score",
    min.cutoff = "q9",
    order = TRUE
  ) &
    NoAxes()

  ggsave(
    filename = file.path(
      output.path,
      "S.Score.png"
    ),
    plot = p,
    width = 7,
    height = 7
  )


  # --------------------------------------------------------------------------
  # G2/M-phase score
  # --------------------------------------------------------------------------

  p <- FeaturePlot(
    data_cluster,
    features = "G2M.Score",
    min.cutoff = "q9",
    order = TRUE
  ) &
    NoAxes()

  ggsave(
    filename = file.path(
      output.path,
      "G2M.Score.png"
    ),
    plot = p,
    width = 7,
    height = 7
  )
}


print(
  "Evaluation of cluster-removal strategies completed."
)
