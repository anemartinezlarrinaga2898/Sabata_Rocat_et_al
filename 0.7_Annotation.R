# SCRIPT: emove some of the clusters as they are noise or not clear what they are 
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 10.07.2024

###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)
library(gridExtra)
library(patchwork)

getPalette <-  colorRampPalette(brewer.pal(8, "Set1"))
col <-  getPalette(20)
path.guardar <- "0.2_SeuratPipeline/Clustering/RemovingCluster"
dir.create(path.guardar)

SeuratPipeline_Subset <- function(data_cluster){
  gene.info.distribution <- summary(Matrix::colSums(data_cluster@assays$RNA@counts[,]>0))
  hvg.number <- round(gene.info.distribution[4]+100)
  data_cluster <- FindVariableFeatures(object = data_cluster,selection.method = "vst", nfeatures = hvg.number)
  data_cluster <- ScaleData(object = data_cluster)
  
  data_cluster <- RunPCA(object = data_cluster)
  
  # Estimate the PCA Dimmensions. 
  pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  dim.final <- min(co1, co2)
  
  data_cluster <- FindNeighbors(object = data_cluster, dims = 1:dim.final)
  data_cluster <- RunUMAP(object = data_cluster, dims = 1:dim.final)
  resolutions <- c(0.1,0.3,0.5,0.7,0.8,0.9)
  data_cluster <- FindClusters(data_cluster, resolution = resolutions)
  return(data_cluster)
}

source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_Annotations_Cluster.R")
source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_SeuratPipeline.R")


Annotation_Function <- function(data,Lista_Markers,path.guardar){
  
  features_df<-do.call(rbind,Lista_Markers)
  cell_type<-names(Lista_Markers)
  
  Lista_Signature <- list()
  
  for(i in seq_along(Lista_Markers)){
    ct<-cell_type[i]
    print(ct)
    m<-features_df$Features[which(features_df$CellType==ct)]
    Lista_Signature[[i]]<-m
    names(Lista_Signature)[i]<-ct
  }
  
  data <- UCell::AddModuleScore_UCell(data, features = Lista_Signature)
  signature_names <- paste0(names(Lista_Signature), "_UCell")
  
  Lista_FP<-list()
  Lista_VL<-list()
  
  for(i in seq_along(signature_names)){
    sig_names<-signature_names[i]
    print(sig_names)
    Fp<-FeaturePlot(data, features = sig_names, min.cutoff = "q9", order = T) &NoAxes()
    Vln<-VlnPlot(data, features =sig_names ,sort = T,cols=col)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    Lista_FP[[i]]<-Fp
    Lista_VL[[i]]<-Vln
  }
  
  x<-do.call("grid.arrange",c(Lista_FP,ncol=5,nrow=2))
  ggsave(plot=x,filename=paste(path.guardar,"FeaturePlot_Annotation.png",sep="/"),width=30,heigh=12)
  
  x<-do.call("grid.arrange",c(Lista_VL,ncol=5,nrow=2))
  ggsave(plot=x,filename=paste(path.guardar,"ViolinPlot_Annotation.png",sep="/"),width=30,heigh=12)
  
}

###############################################################################

data <- readRDS("0.2_SeuratPipeline/Seu.Obj.rds")

dim.final <- Calculo.Dimensiones.PCA(data)
data <- FindNeighbors(object = data, dims = 1:dim.final)
data <- RunUMAP(object = data, dims = 1:dim.final)
resolutions <- c(0.1,0.3,0.5,0.7,0.8,0.9)
data <- FindClusters(data, resolution = resolutions)

data <- SetIdent(data,value="RNA_snn_res.0.8")

EndothelialCells<-data.frame(Features=c("Pecam1","Cdh5"),CellType="EndothelialCells")
Artery<-data.frame(Features=c("Sox17","Dkk2","Efbn2","Sema3g","Fbln5","Hey1","Bmx","Mgp"),CellType="Artery")
Venous<-data.frame(Features=c("Ackr1","Selp","Sele","Icam1","Vwf","Vcam1","Nr2f2","Ephb4"),CellType="Venous")
Capillary<-data.frame(Features=c("Rgcc","Cd36","Kdr","Car2","Cd200","Cd300lg","Giphbp1","Aqp1"),CellType="Capillary")
TipCells<-data.frame(Features=c("Esm1","Dll4","Cxcr4","Col4a2","Col4a1","Sparc","Insr","Angpt2"),CellType="TipCells")
Lymphatics <- data.frame(Features=c("Prox1","Lyve1","Pdpn","Cd63","Flt4","Ccdc3","Mmrn2","Foxp2","Cldn11","Net1","Neo1"),CellType="Lymphatics")
Fenestration <- data.frame(Features=c("Plvap","Plpp3","Igfbp3"),CellType="Fenestration")
Angiogenesis <- data.frame(Features=c("Angpt1","Pdgfc","Vegfc","Vegfb","Lrp1","Tek","Kdr","Angpt2"),CellType="Angiogenesis")
Immature <- data.frame(Features=c("Pdlim1","Rplp0","Gapdh","Igfbp7","Rps2","Rpsa","Rpl12","Aplnr"),CellType="Immature")

Lista_Markers <- list(EndothelialCells=EndothelialCells,Artery=Artery,Venous=Venous,Capillary=Capillary,TipCells=TipCells,Lymphatics=Lymphatics,Fenestration=Fenestration,Angiogenesis=Angiogenesis,Immature=Immature)

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <-  stringr::str_to_title(cc.genes$g2m.genes)

# OPTION 1 ..................................................................... 

clusters <- c("7","9")
data_cluster <- subset(data,idents = clusters,invert=TRUE)
path.guardar <- "0.2_SeuratPipeline/Clustering/RemovingCluster"
path.guardar <- paste(path.guardar,"Remove_Clus_7_9",sep="/")
dir.create(path.guardar)
data_cluster <- SeuratPipeline_Subset(data_cluster)
Annotation_Function(data_cluster,Lista_Markers,path.guardar)
data_cluster <- CellCycleScoring(data_cluster, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(data_cluster,reduction = "umap",group.by = "Phase",label = TRUE, label.size = 5,cols =c("#95BDED","#F5A040","#C8C2F2"),pt.size = 0.5,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phase.png",sep="/"),width = 7,height = 7)

FeaturePlot(data_cluster, features = "S.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score.png",sep="/"),width = 7,height = 7)
FeaturePlot(data_cluster, features = "G2M.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score.png",sep="/"),width = 7,height = 7)

# OPTION 2 ..................................................................... 
clusters <- c("7","9","8")
data_cluster <- subset(data,idents = clusters,invert=TRUE)
path.guardar <- "0.2_SeuratPipeline/Clustering/RemovingCluster"
path.guardar <- paste(path.guardar,"Remove_Clus_7_9_8",sep="/")
dir.create(path.guardar)
data_cluster <- SeuratPipeline_Subset(data_cluster)
Annotation_Function(data_cluster,Lista_Markers,path.guardar)
data_cluster <- CellCycleScoring(data_cluster, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(data_cluster,reduction = "umap",group.by = "Phase",label = TRUE, label.size = 5,cols =c("#95BDED","#F5A040","#C8C2F2"),pt.size = 0.5,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phase.png",sep="/"),width = 7,height = 7)

FeaturePlot(data_cluster, features = "S.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score.png",sep="/"),width = 7,height = 7)
FeaturePlot(data_cluster, features = "G2M.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score.png",sep="/"),width = 7,height = 7)

# OPTION 3 ..................................................................... 
clusters <- c("5","7","9","8")
data_cluster <- subset(data,idents = clusters,invert=TRUE)
path.guardar <- "0.2_SeuratPipeline/Clustering/RemovingCluster"
path.guardar <- paste(path.guardar,"Remove_Clus_7_9_8_5",sep="/")
dir.create(path.guardar)
data_cluster <- SeuratPipeline_Subset(data_cluster)
Annotation_Function(data_cluster,Lista_Markers,path.guardar)
data_cluster <- CellCycleScoring(data_cluster, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(data_cluster,reduction = "umap",group.by = "Phase",label = TRUE, label.size = 5,cols =c("#95BDED","#F5A040","#C8C2F2"),pt.size = 0.5,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phase.png",sep="/"),width = 7,height = 7)

FeaturePlot(data_cluster, features = "S.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score.png",sep="/"),width = 7,height = 7)
FeaturePlot(data_cluster, features = "G2M.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score.png",sep="/"),width = 7,height = 7)





