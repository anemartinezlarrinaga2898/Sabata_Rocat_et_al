# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

################################################################################

# Setting working parameters
source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Harmony",sep="/")
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


