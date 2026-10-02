################################################################################
# PUBLICATION CODE
# FILE: 1.4_Annotate_WT.R
# PURPOSE: Perform preliminary endothelial-state annotation and inspect eGFP/Apln/Esm1 expression.
# NOTE: Original analytical selections and thresholds are retained.
#       Repository-local paths are defined in R/project_config.R.
################################################################################

# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

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
library(plyr)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67/AnnotationsPlots",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-"Harmony_Log_res.0.7"
data@meta.data$Layer_1<- revalue(data@meta.data[[ident]], c("0" = "Cap",
                                                               "1" = "Prolif_1",
                                                               "2" = "Tip",
                                                               "3" = "Prolif_2",
                                                               "4" = "Cap_Artery",
                                                               "5" = "Artery",
                                                               "6" = "Prolif_3",
                                                               "7" = "Angiogenic",
                                                               "8" = "Prolif_4",
                                                               "9" = "Vein",
                                                               "10" = "Remodelling",
                                                               "11" = "Tip"))

levels_Total_layer1 <- c("Vein","Artery","Cap","Cap_Artery","Remodelling","Tip","Angiogenic","Prolif_1","Prolif_2","Prolif_3","Prolif_4")
data$Layer_1<-factor(data$Layer_1,levels=levels_Total_layer1)          

# ----------------------------------------------------
# | Dimplot
# ----------------------------------------------------

ident<-"Layer_1"
DimPlot<-DimPlot(data, group.by = "Layer_1", pt.size = 0.3,label = F,cols=c("#B39DDB","#9FA8DA","#90CAF9","#C5E1A5","#FFF59D","#FFAB91","#F38D30","#FFC685","#916430","#BCAAA4","#CFD8DC"),raster=FALSE)&NoAxes()
file_name<-paste("UMAP_Annot_Clus_",ident,".png",sep="_")
ggsave(plot=DimPlot,filename = paste(path.guardar,file_name,sep="/"),width = 7,height = 7)

saveRDS(data,file.path(path.guardar,"Data_Anotado.rds"))

# ----------------------------------------------------
# | Egfp Levels Definition 
# ----------------------------------------------------

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

DimPlot(data,group.by = "eGFP_Levels",pt.size = 0.5,raster=FALSE,cols=c("blue","grey")) &NoAxes()
ggsave(filename = paste(path.guardar,"UMAP_eGFP_Localization.png",sep="/"),width = 7,height = 7)

matchSCore2::summary_barplot(class.fac = data$eGFP_Levels,obs.fac =data$Layer_1)+scale_fill_manual(values = c("blue","grey"))
ggsave(filename = paste(path.guardar,"BarPlot_eGFP_Distribution.png",sep="/"),width = 6,height = 3)

# ----------------------------------------------------
# | FP Expression
# ----------------------------------------------------

# APLN: 
FeaturePlot(data, features = "Apln", min.cutoff = "q9", order = T)
ggsave(filename = paste(path.guardar,"FP_Apln_Expression.png",sep="/"),width = 7,height = 7)

plot_df <- FetchData(data, vars = c("Apln","eGFP_Levels"))
colnames(plot_df)[1] <- "Expression"
ggplot(plot_df, aes(x = Expression)) +
  # alpha = 0.5 hace que los colores sean semi-transparentes para ver el solapamiento
  geom_density(alpha = 0.5, color = "black", size = 0.3,fill="grey") + 
  # Colores limpios para publicación
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
ggsave(filename = paste(path.guardar,"Density_Apln_Global_Exp.png",sep="/"),width = 6,height = 3)

ggplot(plot_df, aes(x = Expression, fill = eGFP_Levels)) +
  # alpha = 0.5 hace que los colores sean semi-transparentes para ver el solapamiento
  geom_density(alpha = 0.5, color = "black", size = 0.3) + 
  # Colores limpios para publicación
  theme_classic() +
  scale_fill_manual(values = c("blue","grey"))+
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
ggsave(filename = paste(path.guardar,"Density_Apln_Global_Exp_ByEgfP.png",sep="/"),width = 6,height = 3)

# ESM1: 
FeaturePlot(data, features = "Esm1", min.cutoff = "q9", order = T)
ggsave(filename = paste(path.guardar,"FP_Esm1_Expression.png",sep="/"),width = 7,height = 7)

plot_df <- FetchData(data, vars = c("Esm1","eGFP_Levels"))
colnames(plot_df)[1] <- "Expression"
ggplot(plot_df, aes(x = Expression)) +
  # alpha = 0.5 hace que los colores sean semi-transparentes para ver el solapamiento
  geom_density(alpha = 0.5, color = "black", size = 0.3,fill="grey") + 
  # Colores limpios para publicación
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
ggsave(filename = paste(path.guardar,"Density_Esm1_Global_Exp.png",sep="/"),width = 6,height = 3)

ggplot(plot_df, aes(x = Expression, fill = eGFP_Levels)) +
  # alpha = 0.5 hace que los colores sean semi-transparentes para ver el solapamiento
  geom_density(alpha = 0.5, color = "black", size = 0.3) + 
  # Colores limpios para publicación
  theme_classic() +
  scale_fill_manual(values = c("blue","grey"))+
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )
ggsave(filename = paste(path.guardar,"Density_Esm1_Global_Exp_ByEgfP.png",sep="/"),width = 6,height = 3)
