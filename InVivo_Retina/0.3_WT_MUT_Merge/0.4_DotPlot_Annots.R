################################################################################
# SCRIPT: Reference marker DotPlot for combined WT + MUT clustering
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Visualizes a curated marker list from the WT annotation reference file across a user-specified identity in the combined Harmony-integrated dataset.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
#   data/Annot_PCA_24_Info.xlsx (sheet: markers)
#
# OUTPUT:
#   Reference-marker DotPlot in the joint Harmony directory
#
# NOTE:
#   Analytical thresholds, cluster selections, annotations and statistical
#   parameters are preserved from the supplied analysis unless explicitly
#   documented in REVIEW_NOTES.md.
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
data_dir <- file.path(project_root, "data")
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(readxl)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1])) {
  stop("Please provide the metadata identity column as the first command-line argument.")
}
################################################################################

path.obj<-paste(results_root,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"
ident<-args[1]
data<-SetIdent(data,value=ident)
signatures <- read_excel(file.path(data_dir, "Annot_PCA_24_Info.xlsx"),sheet="markers")

# --------------------------------------------------
# | DotPlot
# --------------------------------------------------

DotPlot(data, features = unique(signatures$gene))&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(0.8), face = "plain"), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank())&
  labs(x = "Cell Type", y = "")&
  geom_vline(xintercept = c(5.5,11.5,18.5,25.5), color = "#5D6772", linetype = "dashed")&
  scale_color_gradientn(colours = c("#FFFFFF", "#FEE0B6", "#FC9272", "#DE2D26"))
file_name<-paste("DotPlot_",ident,"_WT_markers.png",sep="")
ggsave(file=paste(path.guardar,file_name,sep="/"),width=12,height=5)
