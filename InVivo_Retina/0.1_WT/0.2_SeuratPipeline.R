################################################################################
# PUBLICATION CODE
# FILE: 0.2_SeuratPipeline.R
# PURPOSE: Run whole-dataset Seurat preprocessing and initial clustering.
# NOTE: Original analytical selections and thresholds are retained.
#       Repository-local paths are defined in R/project_config.R.
################################################################################

# SCRIPT: Normalizacion and Filtering
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

# In this script we are going to filter each of the dataset independtly. We will see if the data has been filtered and
# then depending on that we will decided if the data needs to be filtered or not. However the genes that are detected in less than <10 cells will
# be removed
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
library(writexl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
colors <-  getPalette(200)

library(patchwork)
library(foreach)

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/Total/Seurat",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/Total/Doublets",sep="/")
data <- readRDS(paste(path.obj,"Doublets_Seurat_Total.rds",sep="/"))

# Seurat Pipeline
data <- SeuratPipeline(data)
resolutions <- c(0.1,0.3,0.5,1,1.5,2)
data <- FindClusters(data, resolution = resolutions)

# Save the data
print("Saving data")
saveRDS(data,paste(path.guardar,"SeuratProcess.rds",sep="/")) 
print("Data Saved")

# Plot the UMAPS 
col <-  getPalette(length(unique(data$RNA_snn_res.2)))
umaps <- PlotUMAP(data)
names(umaps) <- paste("Res",resolutions,sep="_")
# Save them as one
umaps.wrap <- wrap_plots(umaps, ncol=4, guides="collect")
ggsave(plot = umaps.wrap,filename = paste(path.guardar,"UMAPS.png",sep="/"),width = 30,height = 20)

# Save each of the plots independently
Guardar <- foreach::foreach(i=1:length(umaps),.final = function(x) setNames(x,names(umaps)))%do%{
  plot.list <- umaps[[i]]
  file.name.saved <- paste(names(umaps)[i],"_UMAP.png",sep="")
  ggsave(plot = plot.list,filename = paste(path.guardar,file.name.saved,sep="/"),width = 10,height = 10)
}

# Plot initial Plots

DimPlot.mt <- FeaturePlot(data,features = "percent.mt",cols = viridis::inferno(25),pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"MT.png",sep="/"),plot=DimPlot.mt,width=10,height=10)
  
DimPlot.counts <-  FeaturePlot(data,features = "nCount_RNA",cols = viridis::inferno(25),pt.size = 1,raster=FALSE)&NoAxes() 
ggsave(filename=paste(path.guardar,"Counts.png",sep="/"),plot=DimPlot.counts,width=10,height=10)

DimPlot.features <-  FeaturePlot(data,features = "nFeature_RNA",cols = viridis::inferno(25),pt.size = 1,raster=FALSE)&NoAxes() 
ggsave(filename=paste(path.guardar,"Features.png",sep="/"),plot=DimPlot.features,width=10,height=10)

col <-  getPalette(length(unique(data$ID)))
DimPlot.ByPatient <- DimPlot(data, group.by = "ID", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.ByPatient.png",sep="/"),plot=DimPlot.ByPatient,width=10,height=10)

DimPlot.Phenotype <- DimPlot(data, group.by = "Phenotype", pt.size = 1,label = FALSE,cols=cols_pheno,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.Phenotype.png",sep="/"),plot=DimPlot.Phenotype,width=10,height=10)

col <-  c("#AADA79","#C2C2C2")
DimPlot.DoubletClassification <- DimPlot(data, group.by = "DoubletClassification_DoubletFinder", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"DimPlot.DoubletClassification.png",sep="/"),plot=DimPlot.DoubletClassification,width=10,height=10)

col <-  getPalette(25)
matchSCore2::summary_barplot(class.fac = data$RNA_snn_res.0.1,obs.fac =data$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res01.png",sep="/"),width = 5,height = 7)

matchSCore2::summary_barplot(class.fac = data$RNA_snn_res.0.3,obs.fac =data$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res03.png",sep="/"),width = 5,height = 7)

FeaturePlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

Fp <- FeaturePlot(data, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes.png",sep="/"),plot=Fp,width=15,height=10)
