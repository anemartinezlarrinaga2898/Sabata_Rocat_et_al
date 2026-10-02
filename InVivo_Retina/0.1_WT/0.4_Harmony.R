################################################################################
# PUBLICATION CODE
# FILE: 0.4_Harmony.R
# PURPOSE: Integrate the complete WT dataset with Harmony using sample ID.
# NOTE: Original analytical selections and thresholds are retained.
#       Repository-local paths are defined in R/project_config.R.
################################################################################

# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 22-11-2023

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
library(harmony)
library(Rcpp)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/Total/Harmony",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")

################################################################################
set.seed(12345)

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/Total/Seurat",sep="/")
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

col <-  c("#AADA79","#C2C2C2")
DimPlot.DoubletClassification <- DimPlot(data.harmony, group.by = "DoubletClassification_DoubletFinder", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.DoubletClassification.png",sep="/"),plot=DimPlot.DoubletClassification,width=10,height=10)

# - ViolinPlot total 
data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.1")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res01.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.3")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res03.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.5")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res05.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="ID")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_ID.png",sep="/"),width = 15,height = 8)
