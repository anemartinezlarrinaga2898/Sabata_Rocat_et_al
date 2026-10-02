################################################################################
# SCRIPT: Integration of DoubletFinder results into the Seurat object
# AUTHOR: Ane Martinez Larrinaga
# DATE: 10-10-2023
#
# DESCRIPTION:
# This script incorporates the previously calculated DoubletFinder results
# into the metadata of each sample-specific Seurat object. Samples are then
# merged into a single Seurat object containing the doublet classification
# information.
#
# INPUT:
#   0.2_SeuratPipeline/Data.Filtered.rds
#   0.2_SeuratPipeline/Doublets/DoubletsRes.rds
#
# OUTPUT:
#   0.2_SeuratPipeline/Data_Doublets.rds
################################################################################


# ------------------------------------------------------------------------------
# 1. Load libraries
# ------------------------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------------------------
# 2. Load DoubletFinder results
# ------------------------------------------------------------------------------

doublet.finder.results <- readRDS(
  "0.2_SeuratPipeline/Doublets/DoubletsRes.rds"
)


# ------------------------------------------------------------------------------
# 3. Load filtered Seurat object
# ------------------------------------------------------------------------------

data <- readRDS(
  "0.2_SeuratPipeline/Data.Filtered.rds"
)


# ------------------------------------------------------------------------------
# 4. Split Seurat object by sample
# ------------------------------------------------------------------------------

data.patient <- SplitObject(
  data,
  split.by = "ID"
)


# ------------------------------------------------------------------------------
# 5. Add DoubletFinder results to sample metadata
# ------------------------------------------------------------------------------

Lista.Patient.Doublets <- vector(
  mode = "list",
  length = length(data.patient)
)

names(Lista.Patient.Doublets) <- names(data.patient)


for (i in seq_along(data.patient)) {

  sample.name <- names(data.patient)[i]

  print(
    paste("Adding DoubletFinder results for:", sample.name)
  )

  # Retrieve sample-specific Seurat object
  sample.data <- data.patient[[i]]

  # Identify corresponding DoubletFinder results
  idx.results.db <- which(
    names(doublet.finder.results) == sample.name
  )

  db.res <- doublet.finder.results[[idx.results.db]]

  # Add DoubletFinder results to Seurat metadata
  sample.data@meta.data <- cbind(
    sample.data@meta.data,
    db.res
  )

  Lista.Patient.Doublets[[i]] <- sample.data
}


# ------------------------------------------------------------------------------
# 6. Merge samples
# ------------------------------------------------------------------------------

Seurat.Object.Total <- merge(
  x = Lista.Patient.Doublets[[1]],
  y = Lista.Patient.Doublets[2:length(Lista.Patient.Doublets)]
)


# ------------------------------------------------------------------------------
# 7. Save Seurat object
# ------------------------------------------------------------------------------

saveRDS(
  Seurat.Object.Total,
  file = "0.2_SeuratPipeline/Data_Doublets.rds"
)

print("DoubletFinder results successfully added to Seurat object.")
