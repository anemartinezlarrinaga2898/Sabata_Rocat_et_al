################################################################################
# SCRIPT: Marker-gene estimation in the eGFP-positive WT + MUT dataset
# AUTHOR: Ane Martinez Larrinaga
# DESCRIPTION:
#   Calculates cluster-wise marker genes for a metadata identity supplied on the
#   command line and generates a DotPlot of the top marker genes.
# INPUT:
#   0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony/Harmony.rds
# ARGUMENT:
#   First trailing argument: metadata column used as the cluster identity.
# OUTPUT:
#   Cluster marker tables (RDS/XLSX) and marker plots.
################################################################################

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), mustWork = FALSE)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- project_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony",sep="/")
dir.create(path.guardar,recursive = TRUE, showWarnings = FALSE)

# Load the script with the functions 
source(file.path(project_root, "utils", "0.3_Util_SeuratPipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide the metadata identity column as the first argument.")
################################################################################

path.obj<-paste(project_root,"0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident <- args[1]
col.idx <- which(colnames(data@meta.data)==ident)
data <- SetIdent(data,value=ident)
clusters<-as.character(unique(data@meta.data[,col.idx]))
clusters<-clusters[order(clusters,decreasing = FALSE)]
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


umap.res.0.5 <- DimPlot(data,reduction = "umap",group.by = "Harmony_Log_res.0.5",label = TRUE, label.size = 5,pt.size = 1,raster=FALSE,split.by="Phenotype")&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.5 ,filename = paste(path.guardar,"Dim_Res05_ByPheno.png",sep="/"),width = 14,height = 7)
