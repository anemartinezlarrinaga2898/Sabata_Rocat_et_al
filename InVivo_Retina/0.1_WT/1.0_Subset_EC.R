################################################################################
# PUBLICATION CODE
# FILE: 1.0_Subset_EC.R
# PURPOSE: Subset candidate endothelial clusters, reprocess them and perform Harmony reintegration.
# NOTE: Original analytical selections and thresholds are retained.
#       Repository-local paths are defined in R/project_config.R.
################################################################################

# SCRIPT: Subset the Endothelial Comparment 
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 30.04.2024

##############################################################################################################################


# Repository configuration ------------------------------------------------------
repo_root <- normalizePath(Sys.getenv("SABATA_REPO_ROOT", unset = "."), mustWork = FALSE)
config_file <- file.path(repo_root, "R", "project_config.R")
if (!file.exists(config_file)) {
  stop("project_config.R not found. Run from the repository root or set SABATA_REPO_ROOT.")
}
source(config_file)

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
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")

# Path_obj
path_obj<-paste(path.guardar_original,"0.4_WT_Analysis/Total/Harmony",sep="/")

##############################################################################################################################

data<-readRDS(paste(path_obj,"ID_LogNormalized.rds",sep="/"))
data<-SetIdent(data,value="Harmony_Log_res.0.3")
# 1) Subset the cells of interest ..............................................................................................

data_interest<-subset(data,idents = c("0","1","2","4","9","11"))
Res<-data.frame(NumberCells=nrow(data_interest@meta.data))
writexl::write_xlsx(Res,paste(path.guardar,"NumberOfCells_BeforeFiltering.xlsx",sep="/"))

# 2) Seurat Pipeline .........................................................................................................
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/Seurat",sep="/")
dir.create(path.guardar, showWarnings = FALSE)

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

# 3) Integration ................................................................................................................................................................................................................................

path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/Harmony",sep="/")
dir.create(path.guardar, showWarnings = FALSE)

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

col <-  getPalette(length(unique(data.harmony$ID)))
DimPlot.ByPatient <- DimPlot(data.harmony, group.by = "ID", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot.ByPatient,width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Phenotype)))
DimPlot.Phenotype.Consensus <- DimPlot(data.harmony, group.by = "Phenotype", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot.Phenotype.Consensus,width=10,height=10)

FeaturePlot(data.harmony, features = c("Pecam1","Cdh5","Vwf","eGFP"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

data.harmony <- SetIdent(data.harmony, value = "Harmony_Log_res.0.3")
VlnPlot(data.harmony, features = c("Pecam1","Cdh5","Vwf","eGFP"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo.png",sep="/"),width = 15,height = 8)
