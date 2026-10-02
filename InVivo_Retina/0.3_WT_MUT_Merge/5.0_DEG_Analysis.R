################################################################################
# SCRIPT: Single-cell differential expression between MUT and WT eGFP-positive endothelial populations
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Performs within-population MUT versus WT differential-expression analysis among eGFP-positive cells using Seurat FindMarkers and exports population-specific results.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/New_Annotations/Data_EC_Annotated.rds
#
# OUTPUT:
#   results/0.6_Join_MUT_WT/DownStream/DEG/New_Annotations/DEG_Mut_WT_eGFP_Pos_Layer_1.{xlsx,rds}
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
library(openxlsx)
library(org.Mm.eg.db)
library(clusterProfiler)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/DownStream/DEG/New_Annotations",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# Load project utility functions 
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
################################################################################

# Data loading and preparation
print("Load data")
# Path obj: 

path.guardar_original <- results_root
path.obj<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/New_Annotations",sep="/")
data <- readRDS(paste(path.obj,"Data_EC_Annotated.rds",sep="/"))
data<-SetIdent(data,value="Layer_1")
print("data loaded")

DefaultAssay(data) <- "RNA"

# Define phenotype levels
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"), labels=c("WT","MUT"))

# Split eGFP SOLO para MUT
data$eGFP_pos <- FetchData(data, vars="eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos, levels=c(FALSE, TRUE), labels=c("FALSE","TRUE"))
data$eGFP_pheno<-paste(data$Phenotype,data$eGFP_pos,sep="_")

ident<-"Layer_1"
data<-SetIdent(data,value=ident)

cluster<-unique(data@active.ident)

Lista_Markers<-list()

for(i in seq_along(cluster)){

    c <- as.character(cluster[i])
    
    # Evitar si el cluster es NA o está vacío
    if(is.na(c) || c == "") {
        next
    }

    data_c <- subset(data, ident = c)
    data_c <- SetIdent(data_c, value = "eGFP_pheno")

    # Contar células por grupo
    cell_counts <- table(Idents(data_c))

    # Comprobar que ambos grupos existen y tienen ≥ 3 células
    if(
        all(c("MUT_TRUE", "WT_TRUE") %in% names(cell_counts)) &&
        cell_counts["MUT_TRUE"] >= 3 &&
        cell_counts["WT_TRUE"] >= 3
    ){

        markers <- FindMarkers(
            data_c,
            ident.1 = "MUT_TRUE",
            ident.2 = "WT_TRUE"
        )

        markers$gene <- rownames(markers)
        markers$cluster <- c
        markers$Classification <- ifelse(markers$avg_log2FC < 0, "Down", "Up")

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


# Comprobar si la lista de marcadores tiene elementos antes de exportar
if(length(Lista_Markers) > 0) {
    file_name <- paste("DEG_Mut_WT_eGFP_Pos_", ident, ".xlsx", sep = "")
    openxlsx::write.xlsx(Lista_Markers, paste(path.guardar, file_name, sep = "/"))
    
    file_name <- paste("DEG_Mut_WT_eGFP_Pos_", ident, ".rds", sep = "")
    saveRDS(Lista_Markers, file.path(path.guardar, file_name))
    print("Results saved successfully")
} else {
    message("No cluster contained enough cells for marker estimation; no output file was generated.")
}
