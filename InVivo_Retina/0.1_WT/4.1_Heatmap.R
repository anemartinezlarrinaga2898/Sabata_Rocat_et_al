################################################################################
# PUBLICATION CODE
# FILE: 5.1_Heatmap.R
# PURPOSE: Generate the final endothelial identity-marker heatmap.
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
library(openxlsx)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))


# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.5_Util_Annotations.R")

Function_Matrix_Exp<-function(data,features_df_filtered){
    clusters<-unique(data@active.ident)
    clusters<-levels(clusters)
    matrix_expression<-matrix(0,ncol=length(clusters),nrow=nrow(features_df_filtered))
    colnames(matrix_expression)<-clusters
    rownames(matrix_expression)<-rownames(features_df_filtered)

    for(i in seq_along(features_df_filtered$Features)){
    f<-features_df_filtered$Features[i]
    print(f)

    idx_f<-which(rownames(matrix_filtered)==f)
    m_f<-matrix_filtered[idx_f,]
    names(m_f)<-colnames(matrix_filtered)

    j<-1

    for(j in seq_along(clusters)){
        c<-clusters[j]
        print(c)

        idx_c<-which(data@active.ident==c)
        cells_id<-rownames(data@meta.data)[idx_c]

        m_f_cells<-m_f[cells_id]
        media<-mean(m_f_cells)
        matrix_expression[i,j]<-media
        }
    }

    return(matrix_expression)
}

args = commandArgs(trailingOnly=TRUE)
################################################################################

# --- Carga y preparación ---
print("Load data")
# Path obj: 

path.guardar_original <- PATHS$results
path.obj<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
data <- readRDS(paste(path.obj,"Data_WT_Anotado.rds",sep="/"))
data<-SetIdent(data,value="AnnotLayer")
print("data loaded")

markers<-readxl::read_excel(PATHS$wt_annotation,sheet="Clean")
colnames(markers)<-c("CellType","Features")

# ----------------------------------------------------
# | Heatmap : 
# ----------------------------------------------------

# ----------------------------------------
# UMAP All populations
cluster_order<-c("Tip Ecs","Angiogenic pre-arterial Ecs","Activated Angiogenic Ecs","Proliferative 1","Proliferative 2","Venous Ecs","Capillary Ecs","Arterial-capillary Ecs","Arterial Ecs")
colors_use<-c("#b37087","#F48FB1","#f4cedb","#b5dcfb","#90caf9","#c5e1a5","#D7B393","#FFE082","#FF8A65")
names(colors_use)<-cluster_order
# ----------------------------------------

features_df<-markers
# 1. Asegurar que es un data.frame clásico (evita el warning de tibble)
features_df <- as.data.frame(features_df)

# 2. Filtrar o limpiar los nombres duplicados (por ejemplo, manteniendo solo la primera aparición de cada gen)
features_df <- features_df[!duplicated(features_df$Features), ]

# 3. Asignar los rownames de forma segura
rownames(features_df) <- features_df$Features

matrix<-data@assays[["RNA"]]@data
genes<-rownames(matrix)
common<-intersect(genes,features_df$Features)
matrix_filtered<-matrix[common,]

features_df_filtered<-features_df[common,]
features_df_filtered<-features_df_filtered[order(features_df_filtered$CellType,decreasing=F),]
rownames(features_df_filtered)<-features_df_filtered$Features

matrix_expression<-Function_Matrix_Exp(data,features_df_filtered)

head(matrix_expression)

matrix_expression <- matrix_expression[, cluster_order]

features_df_filtered_V2 <- data.frame(CellType = features_df_filtered$CellType)
rownames(features_df_filtered_V2) <- rownames(features_df_filtered)

features_df_filtered_V2$CellType <- factor(features_df_filtered_V2$CellType,levels = cluster_order,labels=levels_EC_subset_layer1_8months)
matrix_expression<-t(scale(t(matrix_expression)))

ha <- rowAnnotation(
  CellType = features_df_filtered_V2$CellType,
  col = list(CellType = colors_use),
  show_legend = FALSE
)
col_heatmap <- circlize::colorRamp2(c(-2, 0, 2), c("#1E90FF", "white", "#FF4500"))

pdf(paste(path.guardar,"General_Markers_EC.pdf",sep="/"),width=3,height=8)
Heatmap(matrix_expression,
rect_gp = gpar(col = "white", lwd = 2),
col = col_heatmap,
left_annotation=ha,
cluster_rows = FALSE,
cluster_columns=FALSE,
heatmap_legend_param = list(
        title = "Z-Score"),
column_title = "",
row_names_gp = gpar(fontsize = 8),
column_names_gp = gpar(fontsize = 8),
row_split = features_df_filtered_V2,
row_title = NULL)
dev.off()