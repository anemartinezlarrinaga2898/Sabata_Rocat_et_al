# SCRIPT: Harmony with LogNormalized data
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

################################################################################

# Setting working parameters
source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/FinalPlots",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"
data<-SetIdent(data,value="Harmony_Log_res.0.3")

data$Phenotype<-factor(data$Phenotype,levels=c("WT","MUT"),labels=c("WT","MUT"))

colors_use<-c("#F48FB1","#CE93D8","#9FA8DA","#90CAF9","#80DEEA","#C5E1A5","#FFE082","#FF8A65")

# UMAP Total
DimPlot(data,reduction = "umap",group.by = "Harmony_Log_res.0.3",label = TRUE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types.png",sep="/"),width=7,height=7)
ggsave(paste(path.guardar,"UMAP_All_Cell_Types.pdf",sep="/"),width=7,height=7)

# UMAP Split by Pheno 
DimPlot(data,reduction = "umap",group.by = "Harmony_Log_res.0.3",label = TRUE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE,split.by="Phenotype")&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_SplitByPheno.png",sep="/"),width=14,height=7)
ggsave(paste(path.guardar,"UMAP_SplitByPheno.pdf",sep="/"),width=14,height=7)

# eGFP Levels: 

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 2
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("FALSE","TRUE"),labels=c("FALSE","TRUE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)

DimPlot(data,group.by = "eGFP_pos",pt.size = 0.5,raster=FALSE,cols=c("#E0E0E0","#0D47A1")) &NoAxes()
ggsave(filename=paste(path.guardar,"Dimplot_eGFP.png",sep="/"),width=7,height=7)
ggsave(filename=paste(path.guardar,"Dimplot_eGFP.pdf",sep="/"),width=7,height=7)

DimPlot(data,group.by = "eGFP_pos",pt.size = 0.5,raster=FALSE,split.by="Phenotype",cols=c("#E0E0E0","#0D47A1")) &NoAxes()
ggsave(filename=paste(path.guardar,"Dimplot_eGFP_pos_Split.png",sep="/"),width=14,height=7)
ggsave(filename=paste(path.guardar,"Dimplot_eGFP_pos_Split.pdf",sep="/"),width=14,height=7)

# On the eGFP positive cells the percentage for WT/MUT 

df_clusters <- data@meta.data %>%
  dplyr::mutate(cluster = Idents(data)) %>%
  dplyr::filter(eGFP_pos == "TRUE") %>%   # solo eGFP+
  dplyr::group_by(cluster, Phenotype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  )

ggplot(df_clusters, aes(x = cluster, y = pct, fill = Phenotype)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    x = "Cluster",
    y = "% dentro de eGFP+",
    fill = "Phenotype"
  ) +
  theme_classic()+
  scale_fill_manual(values=c("#E0E0E0","#0D47A1"))
ggsave(filename=paste(path.guardar,"BarPlot_eGFP_Pos_Percentage_Clusters.png",sep="/"),width=5,height=3)
ggsave(filename=paste(path.guardar,"BarPlot_eGFP_Pos_Percentage_Clusters.pdf",sep="/"),width=5,height=3)
