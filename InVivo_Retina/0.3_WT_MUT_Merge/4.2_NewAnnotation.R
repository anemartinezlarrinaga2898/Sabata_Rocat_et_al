################################################################################
# SCRIPT: Revised annotation of the combined WT + MUT endothelial dataset
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Applies the revised endothelial population mapping used for the later differential-expression workflow and saves the updated annotated Seurat object.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
#
# OUTPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/New_Annotations/Data_EC_Annotated.rds and annotation figures
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
library(readxl)
library(plyr)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/New_Annotations",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
################################################################################

path.obj<-paste(results_root,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-"Harmony_Log_res.0.5"
data@meta.data$Layer_1<- revalue(data@meta.data[[ident]], c("0" = "VenousCapillary EC",
                                                               "1" = "Angiogenic EC",
                                                               "2" = "Arterial EC",
                                                               "3" = "Proliferative 1",
                                                               "4" = "Proliferative 2",
                                                               "5" = "Proliferative 1",
                                                               "6" = "Arterial EC",
                                                               "7" = "VenousCapillary EC",
                                                               "8" = "Tip EC"))

data@meta.data$Layer_1 <- as.character(data@meta.data$Layer_1)

idx_artery<-which(data$Harmony_Log_res.0.9=="10")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_artery]<-"Arterial EC"

idx_cap<-which(data$Harmony_Log_res.0.9=="6")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_cap]<-"ArteryCapillary EC"

idx_cap<-which(data$Harmony_Log_res.0.9=="5")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_cap]<-"ArteryCapillary EC"

idx_cap<-which(data$Harmony_Log_res.0.9=="0")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_cap]<-"VenousCapillary EC"

idx_cap<-which(data$Harmony_Log_res.0.9=="9")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_cap]<-"VenousCapillary EC"

data@meta.data$Layer_1[is.na(data@meta.data$Layer_1)] <- "Tip EC"

saveRDS(data,file.path(path.guardar,"Data_EC_Annotated.rds"))

ident<-"Layer_1"
data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"))

# Define colors: 
# ----------------------------------------
# UMAP All populations
cluster_order<-c("Tip EC","Angiogenic EC","Proliferative 1","Proliferative 2","VenousCapillary EC","ArteryCapillary EC","Arterial EC")
colors_use<-c("#b37087","#f4cedb","#b5dcfb","#90caf9","#C5E1A5","#FFE082","#FF8A65")
names(colors_use) <- cluster_order
# ----------------------------------------

data$Layer_1<-factor(data$Layer_1,levels=cluster_order)

# UMAP Total
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = TRUE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_Label.pdf",sep="/"),width=7,height=7)
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = FALSE, cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabel.pdf",sep="/"),width=7,height=7)
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = FALSE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabelNoLegend.pdf",sep="/"),width=7,height=7)

# ----------------------------------------
# Bar Plot of the eGFP proportions 
# ----------------------------------------

data<-SetIdent(data,value="eGFP_pos")
data_eGFP<-subset(data,ident="TRUE")

matchSCore2::summary_barplot(class.fac = data_eGFP$Layer_1,obs.fac =data_eGFP$Phenotype)+scale_fill_manual(values = colors_use)
ggsave(filename=paste(path.guardar,"Barplot_eGFP_WT_MUT.pdf",sep="/"),width=5,height=5)
