# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 22-11-2023

################################################################################

# Setting working parameters
source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(harmony)
library(Rcpp)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.1_SeuratPipeline/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

################################################################################
set.seed(12345)

path.obj<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.1_SeuratPipeline/Seurat",sep="/")
data <- readRDS(paste(path.obj,"SeuratProcess.rds",sep="/"))

DefaultAssay(data)<-"RNA"
integration.var<-"ID"

dim.final <- Calculo.Dimensiones.PCA(data)
print(dim.final)
ElbowPlot(data, ndims = 30, reduction = "pca")
ggsave(paste(path.guardar,"Harmony_ElbowPlot_Log.png",sep="/"),width = 5,height = 5)
    
# Integration with Harmony
data.harmony <- RunHarmony(data,group.by.vars = integration.var,dims = 1:dim.final,plot_convergence=TRUE,reduction.save = "HarmonyLog",assay = "RNA")

data.harmony <- RunUMAP(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final)
data.harmony <- FindNeighbors(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final,graph.name = "Harmony_Log")
resolutions <- c(0.1,0.3,0.5,0.7,0.9,1,1.3,1.5,2)
data.harmony <- FindClusters(data.harmony, resolution = resolutions,graph.name = "Harmony_Log")

# Save the data
print("Savind Data")
file.name<-paste(integration.var,"_LogNormalized.rds",sep="")
saveRDS(data.harmony,paste(path.guardar,file.name,sep="/"))
print("Data Saved")

# 3)Plot the UMAPS and save them
col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.5))+30)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.1))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.1.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.3))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.3", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.3.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.5))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.5", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.5.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.7))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.7", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.7.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.9))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.9", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.9.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.1))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.1.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.1.3))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.1.3", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.1.3.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.1.5))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.1.5", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.1.5.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.2))+30)
DimPlot(data.harmony, group.by = "Harmony_Log_res.2", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.2.png",sep="/"),width=10,height=10)

# 4) Get the umaps with other variable 

col <-  getPalette(length(unique(data.harmony$ID))+30)
DimPlot.ByPatient <- DimPlot(data.harmony, group.by = "ID", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_DimPlot.ByPatient.png",sep="/"),plot=DimPlot.ByPatient,width=20,height=10)

DimPlot.Phenotype <- DimPlot(data.harmony, group.by = "Phenotype", pt.size = 1,label = FALSE,cols=cols_pheno,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.Phenotype.png",sep="/"),plot=DimPlot.Phenotype,width=10,height=10)

col <-  c("#AADA79","#C2C2C2")
DimPlot.DoubletClassification <- DimPlot(data.harmony, group.by = "DoubletClassification_DoubletFinder", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.DoubletClassification.png",sep="/"),plot=DimPlot.DoubletClassification,width=10,height=10)
