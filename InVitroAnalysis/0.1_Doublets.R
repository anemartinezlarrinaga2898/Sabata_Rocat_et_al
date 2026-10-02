################################################################################
# SCRIPT: Doublet detection using DoubletFinder
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# This script performs doublet detection independently for each sample.
# The filtered Seurat object is split by sample ID, each sample is processed
# using the standard Seurat pipeline, and doublets are subsequently estimated
# using DoubletFinder.
#
# INPUT:
#   Data.Filtered.rds
#
# OUTPUT:
#   DoubletsRes.rds
#
# DEPENDENCIES:
#   Util_DoubletDetection.R
#   Util_SeuratPipeline.R
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(DoubletFinder)
library(yaml)


# ------------------------------------------------------------------------------
# 2. Load custom functions
# ------------------------------------------------------------------------------

source("Util_DoubletDetection.R")
source("Util_SeuratPipeline.R")


# ------------------------------------------------------------------------------
# 3. Load filtered Seurat object
# ------------------------------------------------------------------------------

data <- readRDS("Data.Filtered.rds")


# ------------------------------------------------------------------------------
# 4. Create output directory
# ------------------------------------------------------------------------------

path.guardar <- "Doublets"

if (!dir.exists(path.guardar)) {
  dir.create(path.guardar, recursive = TRUE)
}


# ------------------------------------------------------------------------------
# 5. Split Seurat object by sample
# ------------------------------------------------------------------------------

data.patient <- SplitObject(
  data,
  split.by = "ID"
)


# ------------------------------------------------------------------------------
# 6. Process each sample independently
# ------------------------------------------------------------------------------

Patient.Data.Process <- foreach::foreach(
  i = seq_along(data.patient),
  .final = function(x) setNames(x, names(data.patient))
) %do% {

  print(
    paste(
      "Running Seurat pipeline for",
      names(data.patient)[i]
    )
  )

  patient.data <- data.patient[[i]]

  patient.data <- SeuratPipeline(
    patient.data
  )
}


# ------------------------------------------------------------------------------
# 7. Estimate doublets using DoubletFinder
# ------------------------------------------------------------------------------

Doublets.Estimation <- foreach::foreach(
  i = seq_along(Patient.Data.Process),
  .final = function(x) setNames(x, names(Patient.Data.Process))
) %do% {

  print(
    paste(
      "Estimating doublets for",
      names(Patient.Data.Process)[i]
    )
  )

  patient.data <- Patient.Data.Process[[i]]

  patient.data <- DoubletDetection_DF(
    patient.data
  )
}


# ------------------------------------------------------------------------------
# 8. Save results
# ------------------------------------------------------------------------------

saveRDS(
  Doublets.Estimation,
  file = file.path(path.guardar, "DoubletsRes.rds")
)

print("Doublet detection completed.")
