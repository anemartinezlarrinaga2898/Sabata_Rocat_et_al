################################################################################
# PUBLICATION CODE
# FILE: 0.5_Markers.R
# PURPOSE: Estimate cluster markers after whole-dataset Harmony integration.
# NOTE: Original analytical selections and thresholds are retained.
#       Repository-local paths are defined in R/project_config.R.
################################################################################

# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

################################################################################


# Repository configuration ------------------------------------------------------
repo_root <- normalizePath(Sys.getenv("SABATA_REPO_ROOT", unset = "."), mustWork = FALSE)
config_file <- file.path(repo_root, "R", "project_config.R")
if (!file.exists(config_file)) {
  stop("project_config.R not found. Run from the repository root or set SABATA_REPO_ROOT.")
}
source(config_file)

# Setting working parameters
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/Total/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/Total/Harmony",sep="/")
data <- readRDS(paste(path.obj,"ID_LogNormalized.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-args[1]
col.idx <- which(colnames(data@meta.data)==ident)
data <- SetIdent(data,value=ident)
clusters<-as.character(unique(data@meta.data[,col.idx]))
clusters<-clusters[order(clusters,decreasing = F)]
markers<-vector(mode="list",length=length(clusters))
names(markers)<-paste("Clus",clusters,sep="_")

i<-1
for(i in seq_along(clusters)){
    c<-clusters[i]
    idx.c <- data@meta.data[,col.idx][which(data@meta.data[,col.idx]==c)]
    print(c)
    if(length(idx.c)<10){
      print("Not enough cells")
      next
    }else{
      m<-FindMarkers(data,ident.1=c)
      m$gene<-rownames(m)
      m$cluster<-c
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
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)


