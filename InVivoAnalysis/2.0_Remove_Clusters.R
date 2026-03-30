# SCRIPT: Subset the Endothelial Comparment 
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 30.04.2024

##############################################################################################################################

source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(SoupX)
library(harmony)
library(Rcpp)
library(optparse)
library(readxl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# ........................................................................................................
# Stablized the path
# ........................................................................................................

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

# Path_obj
path_obj<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5/Harmony",sep="/")

##############################################################################################################################

data<-readRDS(paste(path_obj,"Harmony.rds",sep="/"))
data<-SetIdent(data,value="Harmony_Log_res.0.3")

# 1) Subset the cells of interest ..............................................................................................

data_interest<-subset(data,idents = c("6"),invert=TRUE)
Res<-data.frame(NumberCells=nrow(data_interest@meta.data))
writexl::write_xlsx(Res,paste(path.guardar,"NumberOfCells_BeforeFiltering.xlsx",sep="/"))

# 2) Seurat Pipeline .........................................................................................................
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Seurat",sep="/")
dir.create(path.guardar,recursive=TRUE)

data_interest <- SeuratPipeline_Subset(data_interest)
resolutions <- c(0.1,0.3,0.5,1,1.5,2)
data_interest <- FindClusters(data_interest, resolution = resolutions)

# Save the data
print("Saving data")
saveRDS(data_interest,paste(path.guardar,"Data_EndothelialCells.rds",sep="/")) 
print("Data Saved")

col <-  getPalette(length(unique(data_interest$RNA_snn_res.2)))
  
umap.res.0.1 <- DimPlot(data_interest,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.1 ,filename = paste(path.guardar,"RNA_Res_01_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.3 <- DimPlot(data_interest,reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.3 ,filename = paste(path.guardar,"RNA_Res_03_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.5 <- DimPlot(data_interest,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.5 ,filename = paste(path.guardar,"RNA_Res_05_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.1 <- DimPlot(data_interest,reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.1 ,filename = paste(path.guardar,"RNA_Res_1_UMAP.png",sep="/"),width = 7,height = 7)

DimPlot <- DimPlot(data_interest,group.by = "ID",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot,width=7,height=7)

DimPlot <- DimPlot(data_interest,group.by = "Phenotype",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot,width=7,height=7)

FeaturePlot(data_interest, features = c("Pecam1","Cdh5","Vwf","eGFP"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

Fp <- FeaturePlot(data_interest, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_interest, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_interest, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_interest, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_interest, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_interest, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data_interest, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes.png",sep="/"),plot=Fp,width=15,height=10)

# 3) Integration ................................................................................................................................................................................................................................

path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE)

integration.var<-"ID"

dim.final <- Calculo.Dimensiones.PCA(data_interest)
print(dim.final)

# Integration with Harmony
data.harmony <- RunHarmony(data_interest,group.by.vars = integration.var,dims = 1:dim.final,plot_convergence=TRUE,reduction.save = "HarmonyLog",assay = "RNA")
data.harmony <- RunUMAP(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final)
data.harmony <- FindNeighbors(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final,graph.name = "Harmony_Log")
resolutions <- c(0.1,0.3,0.5,0.7,0.9,1)
data.harmony <- FindClusters(data.harmony, resolution = resolutions,graph.name = "Harmony_Log")
# Save the data
saveRDS(data.harmony,paste(path.guardar,"Harmony.rds",sep="/"))

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.1))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.1.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.3))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.3", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.3.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.5))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.5", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.5.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.7))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.7", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.7.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.9))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.9", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.9.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.1))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.1.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$ID)))
DimPlot.ByPatient <- DimPlot(data.harmony, group.by = "ID", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot.ByPatient,width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Phenotype)))
DimPlot.Phenotype.Consensus <- DimPlot(data.harmony, group.by = "Phenotype", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot.Phenotype.Consensus,width=10,height=10)

FeaturePlot(data.harmony, features = c("Pecam1","Cdh5","Vwf","eGFP"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

Fp <- FeaturePlot(data.harmony, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data.harmony, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes.png",sep="/"),plot=Fp,width=15,height=10)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.3")
VlnPlot(data.harmony,features=c("Cdh5", "Pecam1", "Cldn5", "Tek","Esm1", "Kdr"),sort = T)
ggsave(filename=paste(path.guardar,"Vln_EC_Res03.png",sep="/"),width=10,height=10)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.3")
VlnPlot(data.harmony,features=c("nFeature_RNA","nCount_RNA","percent.mt"),sort = T)
ggsave(filename=paste(path.guardar,"Vln_EC_Res03_Quality.png",sep="/"),width=15,height=5)
