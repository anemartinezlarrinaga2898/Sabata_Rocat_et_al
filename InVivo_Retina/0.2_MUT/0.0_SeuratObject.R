################################################################################
# SCRIPT: MUT dataset extraction and quality-control assessment
# AUTHOR: Ane Martinez Larrinaga
# DATE: 25-03-2025
#
# DESCRIPTION:
# Extracts the MUT samples from the previously filtered retinal scRNA-seq dataset, evaluates quality-control metrics, reapplies the established QC thresholds, and saves the MUT-specific filtered Seurat object.
#
# INPUT:
#   data/Data.Filtered.rds
#
# OUTPUT:
#   results/0.5_MUT_Analysis/Obj/Data.Filtered.rds
#
# NOTE:
#   Analytical thresholds, clustering resolutions and biological selections
#   are preserved from the original analysis unless explicitly documented in
#   REVIEW_NOTES.md.
################################################################################

# ------------------------------------------------------------------------------
# Project-relative paths
# ------------------------------------------------------------------------------

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
results_root <- file.path(project_root, "results")
utils_dir <- file.path(project_root, "utils")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(harmony)
library(Rcpp)
library(readxl)
library(DropletUtils)
library(DelayedArray)

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.5_MUT_Analysis/Obj",sep="/")
dir.create(path.guardar,recursive=TRUE)

################################################################################

# | List Files 

data_total <- readRDS(file.path(project_root, "data", "Data.Filtered.rds"))
data_total <- SetIdent(data_total, value = "ID")
mut_samples <- c("L1152", "L1155")
data_MUT <- subset(data_total, idents = mut_samples)

# | QC Plotting

# QC violin plots
VlnPlot(data_MUT, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster = FALSE)
ggsave(filename = paste(path.guardar,"QC.png",sep="/"),width = 15,height = 8)

# - Complexity 
Complexity <- data_MUT@meta.data %>% ggplot(aes(x=log10GenesPerUMI, color = Phenotype, fill=Phenotype)) +
    geom_density(alpha = 0.2) +
    theme_classic()
ggsave(filename = paste(path.guardar,"Complexity.png",sep="/"),plot = Complexity,width = 10,height = 10)

# - CellCounts
print("Plotting CellCounts")
CellCounts <- data_MUT@meta.data %>% ggplot(aes(x=ID, fill=ID)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(filename = paste(path.guardar,"CellCounts.png",sep="/"),plot = CellCounts,width = 35,height = 10)

# - UMICount
print("Plotting CellCounts")
UMICount <-data_MUT@meta.data %>%
    ggplot(aes(color=ID, x=nCount_RNA, fill= ID)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
ggsave(filename = paste(path.guardar,"UMICount.png",sep="/"),plot = UMICount,width = 35,height = 10)

# - GenesPerCell
print("Plotting GenesPerCell")
GenesPerCell <-data_MUT@meta.data %>% ggplot(aes(color=ID, x=nFeature_RNA, fill= ID)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 300)
ggsave(filename = paste(path.guardar,"GenesPerCell.png",sep="/"),plot = GenesPerCell,width = 35,height = 10)


# Apply the established QC thresholds
data.filtered <- subset(data_MUT,subset=(nFeature_RNA>=250) & (percent.mt<20) & (log10GenesPerUMI>0.8))
VlnPlot(data.filtered, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster = FALSE)
ggsave(filename = paste(path.guardar,"AfterQC.png",sep="/"),width = 15,height = 8)
saveRDS(data.filtered,paste(path.guardar,"Data.Filtered.rds",sep="/"))
