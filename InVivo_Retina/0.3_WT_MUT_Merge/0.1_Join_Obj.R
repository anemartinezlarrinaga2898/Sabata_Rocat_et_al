################################################################################
# SCRIPT: Merge WT and MUT endothelial datasets and perform joint integration
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Merges the final WT endothelial object with the processed MUT endothelial object, defines eGFP-positive cells, reruns the Seurat workflow, and integrates the combined dataset across sample ID using Harmony.
#
# INPUT:
#   results/0.4_WT_Analysis/.../Final_Dist/Data_WT_Anotado.rds
#   results/0.5_MUT_Analysis/.../Harmony/Harmony.rds
#
# OUTPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Seurat/Data_EndothelialCells.rds
#   results/0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
#
# NOTE:
#   Analytical thresholds, cluster selections, annotations and statistical
#   parameters are preserved from the supplied analysis unless explicitly
#   documented in REVIEW_NOTES.md.
################################################################################

# ------------------------------------------------------------------------------
# Project-relative paths
# ------------------------------------------------------------------------------

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
results_root <- file.path(project_root, "results")
utils_dir <- file.path(project_root, "utils")
data_dir <- file.path(project_root, "data")
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(harmony)
library(Rcpp)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# Load project utility functions 
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
################################################################################

# Load WT Obj

path.obj<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
obj_wt<-readRDS(file.path(path.obj,"Data_WT_Anotado.rds"))

# Load MUT Obj

path.obj<-paste(results_root,"0.5_MUT_Analysis/EC_Subset/RemoveClusters/Remove_Res_05_Clus3679/Harmony",sep="/")
obj_mut <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))

# Merge Obj: 

data <-merge(x = obj_wt,y=obj_mut)

# ------------------------------------------------------------------------------------
# Define eGFP Levels
# ------------------------------------------------------------------------------------

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# ------------------------------------------------------------------------------------
# - Seurat Pipeline 
# ------------------------------------------------------------------------------------

path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Seurat",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

data <- SeuratPipeline_Subset(data)
resolutions <- c(0.1,0.3,0.5,1,1.5,2)
data <- FindClusters(data, resolution = resolutions)

# Save the data
print("Saving data")
saveRDS(data,paste(path.guardar,"Data_EndothelialCells.rds",sep="/")) 
print("Data Saved")

col <-  getPalette(length(unique(data$RNA_snn_res.2)))
  
umap.res.0.1 <- DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.1 ,filename = paste(path.guardar,"RNA_Res_01_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.3 <- DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.3 ,filename = paste(path.guardar,"RNA_Res_03_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.5 <- DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.5 ,filename = paste(path.guardar,"RNA_Res_05_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.1 <- DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.1 ,filename = paste(path.guardar,"RNA_Res_1_UMAP.png",sep="/"),width = 7,height = 7)

DimPlot <- DimPlot(data,group.by = "ID",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot,width=7,height=7)

DimPlot <- DimPlot(data,group.by = "Phenotype",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot,width=7,height=7)

# ------------------------------------------------------------------------------------
# - Harmnoy Pipeline 
# ------------------------------------------------------------------------------------

path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

integration.var<-"ID"

dim.final <- Calculo.Dimensiones.PCA(data)
print(dim.final)

# Integration with Harmony
data.harmony <- RunHarmony(data,group.by.vars = integration.var,dims = 1:dim.final,plot_convergence=TRUE,reduction.save = "HarmonyLog",assay = "RNA")
data.harmony <- RunUMAP(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final)
data.harmony <- FindNeighbors(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final,graph.name = "Harmony_Log")
resolutions <- c(0.1,0.3,0.5,0.7,0.9,1)
data.harmony <- FindClusters(data.harmony, resolution = resolutions,graph.name = "Harmony_Log")
# Save the data
saveRDS(data.harmony,paste(path.guardar,"Harmony.rds",sep="/"))

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.1", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.1.png",sep="/"),width=7,height=7)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.3", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.3.png",sep="/"),width=7,height=7)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.5", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.5.png",sep="/"),width=7,height=7)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.7", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.7.png",sep="/"),width=7,height=7)

DimPlot(data.harmony, group.by = "Harmony_Log_res.0.9", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.0.9.png",sep="/"),width=7,height=7)

DimPlot(data.harmony, group.by = "Harmony_Log_res.1", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Harmony_Log_res.1.png",sep="/"),width=7,height=7)


DimPlot(data.harmony, group.by = "AnnotLayer", pt.size = 1,label = TRUE,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"AnnotLayer_WT.png",sep="/"),width=7,height=7)

# 1. Extraer la columna de metadatos
metadata <- data.harmony@meta.data

# 2. Convertir a factor y asegurar que los NA/NaN sean un nivel
# Reemplazamos los NA por una etiqueta específica o simplemente los agrupamos
metadata$AnnotLayer[is.na(metadata$AnnotLayer)] <- "No_Annot"

# 3. Reordenar los niveles: 
# Ponemos "No_Annot" primero para que se dibuje al principio
niveles <- unique(as.character(metadata$AnnotLayer))
niveles <- c(niveles[niveles != "No_Annot"],"No_Annot")

metadata$AnnotLayer <- factor(metadata$AnnotLayer, levels = niveles)

# 4. Asignar de vuelta al objeto Seurat
data.harmony@meta.data <- metadata

# 5. Ejecutar el plot (ahora los NA están al principio de la lista)
DimPlot(data.harmony, group.by = "AnnotLayer", pt.size = 1, label = TRUE, raster = FALSE) & NoAxes()

# 6. Guardar
ggsave(filename = paste(path.guardar, "AnnotLayer_WT.png", sep = "/"), width = 7, height = 7)

# 1. Crear un subset eliminando los NA
# Usamos subset para quedarnos solo con lo que NO es NA
data.subset <- subset(data.harmony, subset = AnnotLayer != "No_Annot") 
# Nota: Si tus valores son literalmente NA, usa: 
# data.subset <- data.harmony[, !is.na(data.harmony$AnnotLayer)]

# 2. Plotear el subset
DimPlot(data.subset, group.by = "AnnotLayer", pt.size = 1, label = FALSE, raster = FALSE) & NoAxes()

# 3. Guardar
ggsave(filename = paste(path.guardar, "AnnotLayer_WT.png", sep = "/"), width = 7, height = 7)

col <-  getPalette(length(unique(data.harmony$ID)))
DimPlot.ByPatient <- DimPlot(data.harmony, group.by = "ID", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot.ByPatient,width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Phenotype)))
DimPlot.Phenotype.Consensus <- DimPlot(data.harmony, group.by = "Phenotype", pt.size = 1,label = FALSE,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot.Phenotype.Consensus,width=10,height=10)

FeaturePlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"), min.cutoff = "q9", order = TRUE,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

data<-SetIdent(data,value="Harmony_Log_res.0.3")
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo.png",sep="/"),width = 15,height = 8)

data<-SetIdent(data,value="Harmony_Log_res.0.5")
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"),ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo_Res_05.png",sep="/"),width = 15,height = 8)

data<-SetIdent(data,value="Harmony_Log_res.0.7")
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"),ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo_Res_07.png",sep="/"),width = 15,height = 8)

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

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.7")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res07.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="ID")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_ID.png",sep="/"),width = 15,height = 8)

Fp <- FeaturePlot(data.harmony, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data.harmony, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data.harmony, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = TRUE,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes.png",sep="/"),plot=Fp,width=15,height=10)
