# SCRIPT: Mejorar el clustering mediante el aumento de la resolucion y de ver si hay que quitar o no los clusters
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 17-07-2024

###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)

getPalette <-  colorRampPalette(brewer.pal(8, "Set1"))
col <-  getPalette(10)
path.guardar <- "0.2_SeuratPipeline/Clustering"
dir.create(path.guardar)

###############################################################################

data <- readRDS("0.2_SeuratPipeline/Seu.Obj_Remove.rds")

# Increase the resolution ......................................................

resolutions <- c(0.6,0.7,0.8,0.9,1)
data <- FindClusters(data, resolution = resolutions)
saveRDS(data,"0.2_SeuratPipeline/Clustering/Data_HigherClusteringRes.rds")

col <-  getPalette(length(unique(data$RNA_snn_res.1))+2)
DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.6",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.6.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.7",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.7.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.8",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.8.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.9",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.9.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.1.png",sep="/"),width = 10,height = 10)

# Estimate the markers for the different resolutions ...........................

# RES 0.6 ------- ------- ------- ------- ------- ------- ------- ------- ------
ident <- "RNA_snn_res.0.6"
print(ident)
data <- SetIdent(data,value=ident)
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
file_name<-paste("Markers_",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)

# RES 0.7 ------- ------- ------- ------- ------- ------- ------- ------- ------
ident <- "RNA_snn_res.0.7"
print(ident)
data <- SetIdent(data,value=ident)
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
file_name<-paste("Markers_",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)

# RES 0.8 ------- ------- ------- ------- ------- ------- ------- ------- ------
ident <- "RNA_snn_res.0.8"
print(ident)
data <- SetIdent(data,value=ident)
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
file_name<-paste("Markers_",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)

# RES 0.9 ------- ------- ------- ------- ------- ------- ------- ------- ------
ident <- "RNA_snn_res.0.9"
print(ident)
data <- SetIdent(data,value=ident)
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
file_name<-paste("Markers_",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)

# RES 1 ------- ------- ------- ------- ------- ------- ------- ------- ------
ident <- "RNA_snn_res.1"
print(ident)
data <- SetIdent(data,value=ident)
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
file_name<-paste("Markers_",ident,".xlsx",sep="")
openxlsx::write.xlsx(markers,paste(path.guardar,file_name,sep="/"))
markers_plot <- markers %>% arrange(cluster, desc(avg_log2FC)) %>% group_by(cluster) %>% slice_head(n = 5)
features <- unique(markers_plot$gene)
file_name<-paste("DotPlot",ident,".png",sep="")
DotPlot(data, features = features)&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave(file=paste(path.guardar,file_name,sep="/"),width=20,height=5)


