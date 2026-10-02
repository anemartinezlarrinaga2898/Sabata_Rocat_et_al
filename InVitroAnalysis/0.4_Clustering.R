################################################################################
# SCRIPT: Evaluation of higher clustering resolutions
# AUTHOR: Ane Martinez Larrinaga
# DATE: 17-07-2024
#
# DESCRIPTION:
# This script evaluates higher graph-based clustering resolutions in the
# scRNA-seq dataset. Clustering is performed at resolutions ranging from
# 0.6 to 1.0. For each resolution, UMAP representations are generated,
# cluster-specific positive marker genes are identified, and the top five
# markers per cluster are visualized using DotPlots.
#
# INPUT:
#   0.2_SeuratPipeline/Seu.Obj_Remove.rds
#
# OUTPUT:
#   0.2_SeuratPipeline/Clustering/Data_HigherClusteringRes.rds
#   UMAP plots for clustering resolutions 0.6-1.0
#   Marker tables for each clustering resolution
#   DotPlots showing the top five marker genes per cluster
#
# MAIN PARAMETERS:
#   Clustering resolutions: 0.6, 0.7, 0.8, 0.9 and 1.0
#   Marker detection: positive markers only
#   Minimum fraction of expressing cells: 0.25
#   Top markers displayed: 5 genes per cluster
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(openxlsx)


# ------------------------------------------------------------------------------
# 2. Define output directory
# ------------------------------------------------------------------------------

path.guardar <- "0.2_SeuratPipeline/Clustering"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 3. Load Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  "0.2_SeuratPipeline/Seu.Obj_Remove.rds"
)


# ------------------------------------------------------------------------------
# 4. Perform clustering at higher resolutions
# ------------------------------------------------------------------------------

resolutions <- c(
  0.6,
  0.7,
  0.8,
  0.9,
  1.0
)

data <- FindClusters(
  object = data,
  resolution = resolutions
)

saveRDS(
  data,
  file = file.path(
    path.guardar,
    "Data_HigherClusteringRes.rds"
  )
)


# ------------------------------------------------------------------------------
# 5. Visualize clustering resolutions
# ------------------------------------------------------------------------------

getPalette <- colorRampPalette(
  brewer.pal(8, "Set1")
)

cluster.colors <- getPalette(
  length(unique(data$RNA_snn_res.1)) + 2
)


for (res in resolutions) {

  ident <- paste0(
    "RNA_snn_res.",
    format(res, nsmall = 1)
  )

  # Seurat stores resolution 1 as RNA_snn_res.1 rather than RNA_snn_res.1.0
  if (res == 1) {
    ident <- "RNA_snn_res.1"
  }

  print(
    paste("Plotting:", ident)
  )

  p <- DimPlot(
    data,
    reduction = "umap",
    group.by = ident,
    label = TRUE,
    label.size = 5,
    cols = cluster.colors,
    pt.size = 1,
    raster = FALSE
  ) &
    NoAxes()

  # Preserve the original visualization:
  # legend retained for resolution 0.6 and removed for the others
  if (res != 0.6) {
    p <- p & NoLegend()
  }

  ggsave(
    filename = file.path(
      path.guardar,
      paste0(ident, ".png")
    ),
    plot = p,
    width = 10,
    height = 10
  )
}


# ------------------------------------------------------------------------------
# 6. Identify and visualize cluster markers at each resolution
# ------------------------------------------------------------------------------

for (res in resolutions) {

  ident <- paste0(
    "RNA_snn_res.",
    format(res, nsmall = 1)
  )

  if (res == 1) {
    ident <- "RNA_snn_res.1"
  }

  print(
    paste("Identifying markers for:", ident)
  )

  # Set cluster identity
  data <- SetIdent(
    data,
    value = ident
  )


  # --------------------------------------------------------------------------
  # Identify positive cluster markers
  # --------------------------------------------------------------------------

  markers <- FindAllMarkers(
    object = data,
    only.pos = TRUE,
    min.pct = 0.25
  )


  # --------------------------------------------------------------------------
  # Save marker table
  # --------------------------------------------------------------------------

  marker.file <- paste0(
    "Markers_",
    ident,
    ".xlsx"
  )

  openxlsx::write.xlsx(
    markers,
    file.path(
      path.guardar,
      marker.file
    )
  )


  # --------------------------------------------------------------------------
  # Select top five markers per cluster
  # --------------------------------------------------------------------------

  markers_plot <- markers %>%
    arrange(
      cluster,
      desc(avg_log2FC)
    ) %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    ungroup()

  features <- unique(
    markers_plot$gene
  )


  # --------------------------------------------------------------------------
  # Generate marker DotPlot
  # --------------------------------------------------------------------------

  p <- DotPlot(
    data,
    features = features
  ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = rel(1.1),
        face = "plain"
      ),
      axis.text.y = element_text(
        size = rel(1.05),
        face = "plain"
      )
    ) +
    scale_colour_gradientn(
      colours = rev(
        brewer.pal(
          n = 11,
          name = "Spectral"
        )
      )
    )


  # --------------------------------------------------------------------------
  # Save DotPlot
  # --------------------------------------------------------------------------

  dotplot.file <- paste0(
    "DotPlot_",
    ident,
    ".png"
  )

  ggsave(
    filename = file.path(
      path.guardar,
      dotplot.file
    ),
    plot = p,
    width = 20,
    height = 5
  )
}


print("Higher clustering resolution analysis completed.")
