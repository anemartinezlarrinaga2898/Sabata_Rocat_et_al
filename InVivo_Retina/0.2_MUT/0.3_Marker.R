################################################################################
# SCRIPT: Marker identification before Harmony integration
# AUTHOR: Ane Martinez Larrinaga
# DATE: 19-12-2023
#
# DESCRIPTION:
# Identifies cluster markers for a user-specified Seurat identity in the non-integrated MUT dataset and exports per-cluster marker tables together with a top-marker DotPlot.
#
# INPUT:
#   results/0.5_MUT_Analysis/Total/Seurat/SeuratProcess.rds
#
# OUTPUT:
#   Marker tables and DotPlot in results/0.5_MUT_Analysis/Total/Seurat/
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

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.5_MUT_Analysis/Total/Seurat",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load project utility functions
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
################################################################################

path.obj<-paste(results_root,"0.5_MUT_Analysis/Total/Seurat",sep="/")
data <- readRDS(paste(path.obj,"SeuratProcess.rds",sep="/"))
DefaultAssay(data)<-"RNA"

if (length(args) < 1) {
  stop(
    "Missing clustering identity. Run as: Rscript 0.3_Marker.R <metadata_column>"
  )
}

ident <- args[1]

if (!ident %in% colnames(data@meta.data)) {
  stop(paste("Metadata column not found:", ident))
}

col.idx <- which(colnames(data@meta.data)==ident)
data <- SetIdent(data,value=ident)
clusters<-as.character(unique(data@meta.data[,col.idx]))
clusters<-clusters[order(clusters,decreasing = FALSE)]
markers<-vector(mode="list",length=length(clusters))
names(markers)<-paste("Clus",clusters,sep="_")

i<-1
for(i in seq_along(clusters)){
    cluster_id <- clusters[i]
    idx.c <- data@meta.data[,col.idx][which(data@meta.data[,col.idx]==cluster_id)]
    print(cluster_id)
    if(length(idx.c)<10){
      print("Not enough cells")
      next
    }else{
      m<-FindMarkers(data, ident.1 = cluster_id)
      m$gene<-rownames(m)
      m$cluster <- cluster_id
      markers[[i]]<-m
    } 
}   

file_name<-paste("Markers",ident,".rds",sep="")
saveRDS(markers,paste(path.guardar,file_name,sep="/"))
file_name<-paste("Markers",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))

# Save the total markers
markers.totales<-do.call(rbind,markers)
file_name<-paste("Markers_Total",ident,".rds",sep="")
saveRDS(markers.totales,paste(path.guardar,file_name,sep="/"))
file_name<-paste("Markers_Total",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers.totales,paste(path.guardar,file_name,sep="/"))

# DotPlot
markers <- FindAllMarkers(data,only.pos = TRUE,min.pct = 0.25)
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)
