################################################################################
# SCRIPT: Generation and quality control of Seurat objects by sample
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# This script imports Cell Ranger count matrices for individual samples,
# generates a Seurat object for each sample, adds sample metadata and
# quality-control metrics, merges all samples into a single Seurat object,
# and applies cell-level quality-control filtering.
#
# INPUT:
#   0.1_CellRanger/              Cell Ranger count matrices
#   0.0_Summaries/SamplesInfo.xlsx
#
# OUTPUT:
#   0.2_SeuratPipeline/Data.Total.rds
#   0.2_SeuratPipeline/Data.Filtered.rds
#   0.2_SeuratPipeline/QC.png
#   0.2_SeuratPipeline/AfterQC.png
#
# QC THRESHOLDS:
#   nFeature_RNA >= 250
#   percent.mt < 20
#   log10GenesPerUMI > 0.8
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)

# Seurat v4 assay structure was used for this analysis
options(Seurat.object.assay.version = "v4")


# ------------------------------------------------------------------------------
# 2. Define output directory
# ------------------------------------------------------------------------------

path.guardar <- "0.2_SeuratPipeline"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 3. Identify Cell Ranger output directories
# ------------------------------------------------------------------------------

files.to.read <- list.files("0.1_CellRanger")


# ------------------------------------------------------------------------------
# 4. Read Cell Ranger count matrices
# ------------------------------------------------------------------------------

Files <- foreach(
  i = seq_along(files.to.read),
  .final = function(x) setNames(x, files.to.read)
) %do% {

  sample.dir <- files.to.read[i]
  print(paste("Reading:", sample.dir))

  files_path <- file.path(
    "0.1_CellRanger",
    sample.dir
  )

  count_matrix <- Read10X(files_path)
}


# ------------------------------------------------------------------------------
# 5. Load sample metadata
# ------------------------------------------------------------------------------

metadata <- readxl::read_excel(
  "0.0_Summaries/SamplesInfo.xlsx"
)


# ------------------------------------------------------------------------------
# 6. Generate individual Seurat objects and calculate QC metrics
# ------------------------------------------------------------------------------

Seu.QC <- vector(
  mode = "list",
  length = length(Files)
)

names(Seu.QC) <- names(Files)

pattern.mito <- "^mt-"
pattern.ribo <- "^rps"

for (i in seq_along(Files)) {

  cm <- Files[[i]]

  sample.name <- names(Files)[i]
  sample.name <- stringr::str_replace_all(
    sample.name,
    "_output",
    ""
  )

  print(paste("Generating Seurat object:", sample.name))

  metadata.s <- metadata[
    metadata$Sample.Name == sample.name,
  ]

  data <- CreateSeuratObject(
    counts = cm,
    project = "mgraupera_10"
  )

  # Add sample information
  data$ID <- metadata.s$Sample.Name
  data$Phenotype <- metadata.s$Phenotype

  # Calculate QC metrics
  data$log10GenesPerUMI <-
    log10(data$nFeature_RNA) / log10(data$nCount_RNA)

  data$percent.mt <- PercentageFeatureSet(
    data,
    pattern = pattern.mito
  )

  data$percent.mt.div100 <- data$percent.mt / 100

  Seu.QC[[i]] <- data
}


# ------------------------------------------------------------------------------
# 7. Merge samples
# ------------------------------------------------------------------------------

data.total <- merge(
  x = Seu.QC[[1]],
  y = Seu.QC[2:length(Seu.QC)],
  add.cell.ids = names(Files)
)

data.total <- SetIdent(
  data.total,
  value = "ID"
)

saveRDS(
  data.total,
  file = file.path(path.guardar, "Data.Total.rds")
)


# ------------------------------------------------------------------------------
# 8. Visualize QC metrics before filtering
# ------------------------------------------------------------------------------

qc.before <- VlnPlot(
  data.total,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  log = TRUE,
  raster = FALSE
)

ggsave(
  filename = file.path(path.guardar, "QC.png"),
  plot = qc.before,
  width = 15,
  height = 8
)


# ------------------------------------------------------------------------------
# 9. Apply cell-level quality-control filtering
# ------------------------------------------------------------------------------

data.filtered <- subset(
  data.total,
  subset =
    nFeature_RNA >= 250 &
    percent.mt < 20 &
    log10GenesPerUMI > 0.8
)


# ------------------------------------------------------------------------------
# 10. Visualize QC metrics after filtering
# ------------------------------------------------------------------------------

qc.after <- VlnPlot(
  data.filtered,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3,
  log = TRUE,
  raster = FALSE
)

ggsave(
  filename = file.path(path.guardar, "AfterQC.png"),
  plot = qc.after,
  width = 15,
  height = 8
)


# ------------------------------------------------------------------------------
# 11. Save filtered Seurat object
# ------------------------------------------------------------------------------

saveRDS(
  data.filtered,
  file = file.path(path.guardar, "Data.Filtered.rds")
)

print("Seurat object generation and quality-control filtering completed.")
