################################################################################
# SCRIPT: eGFP-positive endothelial cell subset, reclustering and Harmony integration
# AUTHOR: Ane Martinez Larrinaga
# DESCRIPTION:
#   Selects eGFP-positive endothelial cells from the combined WT + MUT dataset,
#   performs Seurat processing and clustering, integrates samples with Harmony,
#   and generates QC and marker-expression visualizations.
# INPUT:
#   0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
# OUTPUT:
#   Seurat and Harmony objects plus diagnostic plots under
#   0.6_Join_MUT_WT/EC_Subset/eGFP/.
# IMPORTANT:
#   The eGFP-positive threshold is preserved exactly as in the original analysis
#   (eGFP > 1.5). See REVIEW_NOTES.md for historical-code issues that were not
#   silently changed.
################################################################################

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), mustWork = FALSE)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(plyr)
library(harmony)
library(Rcpp)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- project_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/eGFP",sep="/")
dir.create(path.guardar,recursive = TRUE, showWarnings = FALSE)
source(file.path(project_root, "utils", "0.3_Util_SeuratPipeline.R"))

################################################################################

path.obj<-paste(project_root,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

data<-SetIdent(data,value="eGFP_Levels")
data_subset<-subset(data,idents = c("eGFP +"))


# ------------------------------------------------------------------------------------
# - Seurat Pipeline 
# ------------------------------------------------------------------------------------

path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/eGFP/Seurat",sep="/")
dir.create(path.guardar,recursive = TRUE, showWarnings = FALSE)

data <- SeuratPipeline_Subset(data_subset)
resolutions <- c(0.1,0.3,0.5,1,1.5,2)
data_subset <- FindClusters(data_subset, resolution = resolutions)

# Save the data
print("Saving data")
saveRDS(data_subset,paste(path.guardar,"Data_EndothelialCells.rds",sep="/")) 
print("Data Saved")

col <-  getPalette(length(unique(data_subset$RNA_snn_res.2)))
  
umap.res.0.1 <- DimPlot(data_subset,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.1 ,filename = paste(path.guardar,"RNA_Res_01_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.3 <- DimPlot(data_subset,reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.3 ,filename = paste(path.guardar,"RNA_Res_03_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.0.5 <- DimPlot(data_subset,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.0.5 ,filename = paste(path.guardar,"RNA_Res_05_UMAP.png",sep="/"),width = 7,height = 7)

umap.res.1 <- DimPlot(data_subset,reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(plot = umap.res.1 ,filename = paste(path.guardar,"RNA_Res_1_UMAP.png",sep="/"),width = 7,height = 7)

DimPlot <- DimPlot(data_subset,group.by = "ID",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"ID.png",sep="/"),plot=DimPlot,width=7,height=7)

DimPlot <- DimPlot(data_subset,group.by = "Phenotype",pt.size = 1,raster=FALSE) &NoAxes()
ggsave(filename=paste(path.guardar,"Phenotype.png",sep="/"),plot=DimPlot,width=7,height=7)

# ------------------------------------------------------------------------------------
# - Harmnoy Pipeline 
# ------------------------------------------------------------------------------------

path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony",sep="/")
dir.create(path.guardar,recursive = TRUE, showWarnings = FALSE)

integration.var<-"ID"

dim.final <- Calculo.Dimensiones.PCA(data_subset)
print(dim.final)

# Integration with Harmony
data.harmony <- RunHarmony(data_subset,group.by.vars = integration.var,dims = 1:dim.final,plot_convergence=TRUE,reduction.save = "HarmonyLog",assay = "RNA")
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
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"), ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo.png",sep="/"),width = 15,height = 8)

data<-SetIdent(data,value="Harmony_Log_res.0.5")
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"),ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo_Res_05.png",sep="/"),width = 15,height = 8)

data<-SetIdent(data,value="Harmony_Log_res.0.7")
VlnPlot(data, features = c("Pecam1","Cdh5","Vwf","eGFP"),ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"Vln_Markers_Endo_Res_07.png",sep="/"),width = 15,height = 8)

# - ViolinPlot total 
data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.1")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res01.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.3")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res03.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.5")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res05.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="Harmony_Log_res.0.7")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC_Res07.png",sep="/"),width = 15,height = 8)

data.harmony<-SetIdent(data.harmony,value="ID")
VlnPlot(data.harmony, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = TRUE,raster=FALSE)
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


