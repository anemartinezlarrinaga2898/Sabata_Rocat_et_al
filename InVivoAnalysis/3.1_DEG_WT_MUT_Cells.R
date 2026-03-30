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
library(openxlsx)
library(org.Mm.eg.db)
library(clusterProfiler)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03/Downstream/DEG/Mut_WT",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Subset_Clus_Res_03/Harmony",sep="/")
data <- readRDS(paste(path.data,"Harmony.rds",sep="/"))

DefaultAssay(data)<-"RNA"
data$Phenotype<-factor(data$Phenotype,levels=c("WT","MUT"),labels=c("WT","MUT"))

ident<-args[1]
data<-SetIdent(data,value=ident)
cluster<-unique(data@active.ident)

Lista_Markers<-list()
Lista_GOTerms<-list()

for (i in seq_along(cluster)) {

  c <- cluster[i]

  data_c <- subset(data, ident = c)
  data_c <- SetIdent(data_c, value = "Phenotype")

  # Contar células por grupo
  cell_counts <- table(Idents(data_c))

  # Chequeo mínimo de tamaño
  if (
    all(c("MUT", "WT") %in% names(cell_counts)) &&
    cell_counts["MUT"] >= 3 &&
    cell_counts["WT"] >= 3
  ) {

    markers <- FindMarkers(
      data_c,
      ident.1 = "MUT",
      ident.2 = "WT"
    )

    markers$gene <- rownames(markers)
    markers$cluster <- c
    markers$Classification <- ifelse(
      markers$avg_log2FC < 0,
      "Down",
      "Up"
    )

    Lista_Markers[[length(Lista_Markers) + 1]] <- markers
    names(Lista_Markers)[length(Lista_Markers)] <- paste("Cluster", c, sep = "_")

  } else {

    message(
      paste(
        "Skipping cluster", c,
        "- insufficient cells:",
        paste(names(cell_counts), cell_counts, collapse = ", ")
      )
    )
  }
}

file_name<-paste("DEG_Mut_WT_",ident,".xlsx",sep="")
openxlsx::write.xlsx(Lista_Markers,paste(path.guardar,file_name,sep="/"))
