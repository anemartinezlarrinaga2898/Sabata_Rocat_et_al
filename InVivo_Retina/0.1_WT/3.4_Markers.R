################################################################################
# PUBLICATION CODE
# FILE: 3.4_Markers.R
# PURPOSE: Estimate markers after the sample-specific exclusion.
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
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
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



