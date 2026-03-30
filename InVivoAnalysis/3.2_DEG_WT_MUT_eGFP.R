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
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Downstream/DEG/eGFP_onlyMut_wholeWT",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

# --- Carga y preparación ---
path.data <- paste(Lista_Paths_Main$Path_Guardar_Retinas,
                   "0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony",
                   sep="/")
data <- readRDS(paste(path.data, "Harmony.rds", sep="/"))

DefaultAssay(data) <- "RNA"

# Asegurar niveles
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"), labels=c("WT","MUT"))

# Split eGFP SOLO para MUT
data$eGFP_pos <- FetchData(data, vars="eGFP")$eGFP > 2
data$eGFP_pos <- factor(data$eGFP_pos, levels=c(FALSE, TRUE), labels=c("FALSE","TRUE"))

# Grupo final:
# - WT => "WT_ALL"
# - MUT => "MUT_TRUE" o "MUT_FALSE"
data$group_compare <- ifelse(
  data$Phenotype == "WT",
  "WT_ALL",
  paste0("MUT_", as.character(data$eGFP_pos))
)
data$group_compare <- factor(data$group_compare, levels=c("WT_ALL","MUT_FALSE","MUT_TRUE"))

# --- Identidad base ---
ident <- args[1]
data <- SetIdent(data, value = ident)
cluster <- unique(data@active.ident)

Lista_Markers <- list()

for(i in seq_along(cluster)){

  c <- cluster[i]

  data_c <- subset(data, ident = c)
  data_c <- SetIdent(data_c, value = "group_compare")

  cell_counts <- table(Idents(data_c))

  # Comparación 1: MUT_TRUE vs WT_ALL
  if(all(c("MUT_TRUE","WT_ALL") %in% names(cell_counts)) &&
     cell_counts["MUT_TRUE"] >= 3 &&
     cell_counts["WT_ALL"] >= 3){

    markers_true <- FindMarkers(
      data_c,
      ident.1 = "MUT_TRUE",
      ident.2 = "WT_ALL"
    )

    markers_true$gene <- rownames(markers_true)
    markers_true$cluster <- c
    markers_true$comparison <- "MUT_TRUE_vs_WT_ALL"
    markers_true$Classification <- ifelse(markers_true$avg_log2FC < 0, "Down", "Up")

    Lista_Markers[[length(Lista_Markers) + 1]] <- markers_true
    names(Lista_Markers)[length(Lista_Markers)] <- paste("Cluster", c, "MUT_TRUE_vs_WT_ALL", sep = "_")

  } else {
    message(paste("Skipping cluster", c, "for MUT_TRUE vs WT_ALL - counts:",
                  paste(names(cell_counts), cell_counts, collapse=", ")))
  }

  # Comparación 2: MUT_FALSE vs WT_ALL
  if(all(c("MUT_FALSE","WT_ALL") %in% names(cell_counts)) &&
     cell_counts["MUT_FALSE"] >= 3 &&
     cell_counts["WT_ALL"] >= 3){

    markers_false <- FindMarkers(
      data_c,
      ident.1 = "MUT_FALSE",
      ident.2 = "WT_ALL"
    )

    markers_false$gene <- rownames(markers_false)
    markers_false$cluster <- c
    markers_false$comparison <- "MUT_FALSE_vs_WT_ALL"
    markers_false$Classification <- ifelse(markers_false$avg_log2FC < 0, "Down", "Up")

    Lista_Markers[[length(Lista_Markers) + 1]] <- markers_false
    names(Lista_Markers)[length(Lista_Markers)] <- paste("Cluster", c, "MUT_FALSE_vs_WT_ALL", sep = "_")

  } else {
    message(paste("Skipping cluster", c, "for MUT_FALSE vs WT_ALL - counts:",
                  paste(names(cell_counts), cell_counts, collapse=", ")))
  }
}

file_name <- paste("DEG_MUT_eGFPsplit_vs_WTALL_", ident, ".xlsx", sep="")
openxlsx::write.xlsx(Lista_Markers, paste(path.guardar, file_name, sep="/"))
