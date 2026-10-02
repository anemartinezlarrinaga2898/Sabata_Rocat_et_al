################################################################################
# SCRIPT: UCell signature scoring and cell-cycle assessment
# AUTHOR: Ane Martinez Larrinaga
# DATE: 29-07-2024
#
# DESCRIPTION:
# This script evaluates predefined gene signatures across the endothelial
# scRNA-seq dataset after removal of non-endothelial/low-quality clusters.
#
# Gene signatures are imported from Signatures.xlsx and scored at single-cell
# level using UCell. Signature scores are visualized across the UMAP and
# clustering structure using FeaturePlots and violin plots. Expression of
# individual genes belonging to each signature is also visualized.
#
# In addition, cell-cycle scores are calculated using Seurat CellCycleScoring
# to evaluate S and G2/M transcriptional programs and their distribution
# across experimental phenotypes.
#
# INPUT:
#   0.2_SeuratPipeline/Seu.Obj_Remove.rds
#   0.0_Summaries/Signatures.xlsx
#
# OUTPUT:
#   UCell signature FeaturePlots
#   UCell signature violin plots
#   Individual signature-gene FeaturePlots
#   Cell-cycle phase UMAP
#   S-phase and G2/M score FeaturePlots
#
# CLUSTER IDENTITY:
#   RNA_snn_res.0.3
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(gridExtra)
library(UCell)


# ------------------------------------------------------------------------------
# 2. Define paths and plotting parameters
# ------------------------------------------------------------------------------

path.guardar <- "0.2_SeuratPipeline/Annotations"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}

getPalette <- colorRampPalette(
  brewer.pal(8, "Set1")
)

col <- getPalette(20)


# ------------------------------------------------------------------------------
# 3. Load endothelial Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  "0.2_SeuratPipeline/Seu.Obj_Remove.rds"
)

data$Phenotype <- factor(
  data$Phenotype,
  levels = c("OB_et", "OB_40HT"),
  labels = c("OB_et", "OB_40HT")
)

data <- SetIdent(
  data,
  value = "RNA_snn_res.0.3"
)


# ------------------------------------------------------------------------------
# 4. Load predefined gene signatures
# ------------------------------------------------------------------------------

Signatures <- read_excel(
  "0.0_Summaries/Signatures.xlsx"
)

Signatures_List <- split(
  Signatures,
  Signatures$Signature
)

names(Signatures_List) <- stringr::str_replace_all(
  names(Signatures_List),
  " ",
  "_"
)

names(Signatures_List) <- stringr::str_to_title(
  names(Signatures_List)
)


# ------------------------------------------------------------------------------
# 5. Generate gene-signature list
# ------------------------------------------------------------------------------

Lista_Signature <- list()

for (i in seq_along(Signatures_List)) {

  signature.data <- Signatures_List[[i]]

  genes <- stringr::str_to_title(
    signature.data$Gene
  )

  Lista_Signature[[i]] <- genes

  names(Lista_Signature)[i] <-
    names(Signatures_List)[i]
}


# ------------------------------------------------------------------------------
# 6. Calculate UCell signature scores
# ------------------------------------------------------------------------------

data <- UCell::AddModuleScore_UCell(
  data,
  features = Lista_Signature
)

signature_names <- paste0(
  names(Lista_Signature),
  "_UCell"
)


# ------------------------------------------------------------------------------
# 7. Generate combined signature plots
# ------------------------------------------------------------------------------

Lista_FP <- list()
Lista_VL <- list()

for (i in seq_along(signature_names)) {

  signature <- signature_names[i]

  print(
    paste("Plotting signature:", signature)
  )

  Lista_FP[[i]] <- FeaturePlot(
    data,
    features = signature,
    min.cutoff = "q9",
    order = TRUE
  ) &
    NoAxes()

  Lista_VL[[i]] <- VlnPlot(
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


# Combined FeaturePlots
p <- do.call(
  "grid.arrange",
  c(
    Lista_FP,
    ncol = 5,
    nrow = 2
  )
)

ggsave(
  filename = file.path(
    path.guardar,
    "FeaturePlot_Annotation.png"
  ),
  plot = p,
  width = 30,
  height = 12
)


# Combined violin plots
p <- do.call(
  "grid.arrange",
  c(
    Lista_VL,
    ncol = 5,
    nrow = 2
  )
)

ggsave(
  filename = file.path(
    path.guardar,
    "ViolinPlot_Annotation.png"
  ),
  plot = p,
  width = 30,
  height = 12
)


# ------------------------------------------------------------------------------
# 8. Generate individual plots for each UCell signature
# ------------------------------------------------------------------------------

for (signature in signature_names) {

  print(
    paste("Saving individual signature:", signature)
  )


  # UMAP FeaturePlot
  p <- FeaturePlot(
    data,
    features = signature,
    min.cutoff = "q9",
    order = TRUE
  ) &
    NoAxes()

  ggsave(
    filename = file.path(
      path.guardar,
      paste0("Fp_", signature, ".png")
    ),
    plot = p,
    width = 5,
    height = 5
  )


  # Violin plot
  p <- VlnPlot(
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

  ggsave(
    filename = file.path(
      path.guardar,
      paste0("Vln_", signature, ".png")
    ),
    plot = p,
    width = 5,
    height = 5
  )
}


# ------------------------------------------------------------------------------
# 9. Visualize individual genes from each signature
# ------------------------------------------------------------------------------

for (i in seq_along(Lista_Signature)) {

  signature.genes <- Lista_Signature[[i]]
  signature.name <- names(Lista_Signature)[i]

  print(
    paste("Individual genes for:", signature.name)
  )

  signature.path <- file.path(
    path.guardar,
    signature.name
  )

  if (!dir.exists(signature.path)) {
    dir.create(
      signature.path,
      recursive = TRUE
    )
  }


  for (gene in signature.genes) {

    # Only plot genes detected in the RNA assay
    if (gene %in% rownames(data[["RNA"]])) {

      p <- FeaturePlot(
        data,
        features = gene,
        min.cutoff = "q9",
        order = TRUE
      ) &
        NoAxes()

      ggsave(
        filename = file.path(
          signature.path,
          paste0(gene, ".png")
        ),
        plot = p,
        width = 5,
        height = 5
      )
    }
  }
}


# ------------------------------------------------------------------------------
# 10. Cell-cycle scoring
# ------------------------------------------------------------------------------

s.genes <- stringr::str_to_title(
  cc.genes$s.genes
)

g2m.genes <- stringr::str_to_title(
  cc.genes$g2m.genes
)

data <- CellCycleScoring(
  data,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)


# ------------------------------------------------------------------------------
# 11. Visualize phenotype and cell-cycle phase
# ------------------------------------------------------------------------------

phase.cols <- c(
  "#95BDED",
  "#F5A040",
  "#C8C2F2"
)


# Phenotype split by cell-cycle phase
p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "Phenotype",
  split.by = "Phase",
  label = FALSE,
  cols = scales::alpha(phase.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "Phenotype_SplitByPhase.png"
  ),
  plot = p,
  width = 15,
  height = 7
)


# Phenotype distribution
p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "Phenotype",
  label = FALSE,
  cols = scales::alpha(phase.cols, 0.66),
  pt.size = 1,
  raster = FALSE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "Phenotype.png"
  ),
  plot = p,
  width = 7,
  height = 7
)


# Cell-cycle phase
p <- DimPlot(
  data,
  reduction = "umap",
  group.by = "Phase",
  label = TRUE,
  label.size = 5,
  cols = phase.cols,
  pt.size = 0.5,
  raster = FALSE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "Phase.png"
  ),
  plot = p,
  width = 7,
  height = 7
)


# ------------------------------------------------------------------------------
# 12. Visualize S and G2/M scores
# ------------------------------------------------------------------------------

p <- FeaturePlot(
  data,
  features = "S.Score",
  min.cutoff = "q9",
  order = TRUE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "S.Score.png"
  ),
  plot = p,
  width = 7,
  height = 7
)


p <- FeaturePlot(
  data,
  features = "G2M.Score",
  min.cutoff = "q9",
  order = TRUE
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "G2M.Score.png"
  ),
  plot = p,
  width = 7,
  height = 7
)


# ------------------------------------------------------------------------------
# 13. Visualize cell-cycle scores by phenotype
# ------------------------------------------------------------------------------

p <- FeaturePlot(
  data,
  features = "S.Score",
  min.cutoff = "q9",
  order = TRUE,
  split.by = "Phenotype"
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "S.Score_SplitByPhenotype.png"
  ),
  plot = p,
  width = 10,
  height = 7
)


p <- FeaturePlot(
  data,
  features = "G2M.Score",
  min.cutoff = "q9",
  order = TRUE,
  split.by = "Phenotype"
) &
  NoAxes()

ggsave(
  filename = file.path(
    path.guardar,
    "G2M.Score_SplitByPhenotype.png"
  ),
  plot = p,
  width = 10,
  height = 7
)


print(
  "UCell signature scoring and cell-cycle assessment completed."
)
