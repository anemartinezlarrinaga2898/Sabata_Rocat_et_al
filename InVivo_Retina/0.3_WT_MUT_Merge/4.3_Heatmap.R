################################################################################
# SCRIPT: Expression similarity heatmaps across WT and MUT endothelial populations
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Calculates average-expression profiles by phenotype and endothelial population and visualizes correlation-based hierarchical clustering for all cells and the eGFP-positive subset.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Annotations/Data_EC_Annotated.rds
#
# OUTPUT:
#   Average-expression hierarchical-clustering heatmaps
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
library(ComplexHeatmap)
library(circlize)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/DownStream/Similarity",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)
################################################################################

path.guardar_original <- results_root
path.obj<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Annotations",sep="/")
data <- readRDS(paste(path.obj,"Data_EC_Annotated.rds",sep="/"))
data<-SetIdent(data,value="Layer_1")
print("data loaded")

data$Condition_Cluster<-paste(data@active.ident,data$Phenotype,sep="_")
ident<-"Condition_Cluster"

# --------------------------------------------------------------------------------------------------------------------------------------------
# Heatmap: WT + MUT together with LABELS 
# --------------------------------------------------------------------------------------------------------------------------------------------

# --- PARÁMETROS QUE PUEDES EDITAR ---
file_name<-paste(ident,"HC_AverageExpression_Total.pdf",sep="_")
outfile <- file.path(path.guardar, file_name)
group_by <- ident            # columna de meta que define columnas del heatmap
palette_cont <- "RdBu"                 # "viridis" | "RdBu" | "PuOr"
cap_limits <- c(-2.5, 2.5)             # capping del z-score
top_var_rows <- 1500                   # opcional: top genes por varianza (mejora señal/ruido)
# ------------------------------------

suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize); library(RColorBrewer)
})

# 1) matriz de expresión media (genes x grupos)
hvg <- VariableFeatures(data)
avg <- AverageExpression(data, features = hvg, group.by = group_by, assays = "RNA", slot = "data")$RNA
avg <- avg[apply(avg, 1, var) > 0, , drop = FALSE]

# 2) seleccionar top-N por varianza (opcional)
if (!is.null(top_var_rows) && top_var_rows < nrow(avg)) {
  v <- apply(avg, 1, var)
  avg <- avg[order(v, decreasing = TRUE)[seq_len(top_var_rows)], , drop = FALSE]
}

# 3) z-score por fila y capping
z <- t(scale(t(as.matrix(avg))))
z[z < cap_limits[1]] <- cap_limits[1]
z[z > cap_limits[2]] <- cap_limits[2]

# 4) función de color
make_col_fun <- function(low, high, palette=c("viridis","RdBu","PuOr")[2]) {
  if (palette=="viridis") cols <- viridisLite::viridis(11) else
  if (palette=="RdBu")   cols <- rev(brewer.pal(11,"RdBu")) else
                           cols <- colorRampPalette(brewer.pal(9,"PuOr"))(11)
  circlize::colorRamp2(seq(low, high, length.out=length(cols)), cols)
}
col_fun <- make_col_fun(cap_limits[1], cap_limits[2], palette_cont)

# 5) distancia tipo 1 - correlación (mejor para perfiles de expresión)
dist_cols <- as.dist(1 - cor(z, method="pearson"))
dist_rows <- as.dist(1 - cor(t(z), method="pearson"))

ht <- Heatmap(
  z, name="Z-score", col=col_fun,
  clustering_distance_columns = dist_cols,
  clustering_distance_rows    = dist_rows,
  clustering_method_columns   = "average",
  clustering_method_rows      = "average",
  show_row_names = FALSE, show_column_names = TRUE,
  column_title = "Average expression (z-score per gene)",
  border=TRUE, rect_gp = gpar(col="#EEEEEE", lwd=0.4),
  heatmap_legend_param = list(at=seq(cap_limits[1],cap_limits[2],by=1),
                              title_gp=gpar(fontface="bold"), labels_gp=gpar(fontsize=8)),
  use_raster=TRUE, raster_quality=3
)

# 6) exportación vectorial
dims_w <- 3.5 + 0.22 * ncol(z); dims_h <- 4.0 + 0.18 * nrow(z)
dims_w <- max(6, min(dims_w, 16)); dims_h <- max(6, min(dims_h, 20))

grDevices::cairo_pdf(outfile, width=5, height=10, family="Arial")
draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
message("Guardado: ", outfile)

# --------------------------------------------------------------------------------------------------------------------------------------------
# Heatmap: WT + MUT together with LABELS only in EGFP+
# --------------------------------------------------------------------------------------------------------------------------------------------

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"))

data<-SetIdent(data,value="eGFP_pos")
data_eGFP<-subset(data,ident="TRUE")

ident<-"Condition_Cluster"

# --- PARÁMETROS QUE PUEDES EDITAR ---
file_name<-paste(ident,"HC_AverageExpression_eGFP.pdf",sep="_")
outfile <- file.path(path.guardar, file_name)
group_by <- ident            # columna de meta que define columnas del heatmap
palette_cont <- "RdBu"                 # "viridis" | "RdBu" | "PuOr"
cap_limits <- c(-2.5, 2.5)             # capping del z-score
top_var_rows <- 1500                   # opcional: top genes por varianza (mejora señal/ruido)
# ------------------------------------

suppressPackageStartupMessages({
})

# 1) matriz de expresión media (genes x grupos)
hvg <- VariableFeatures(data_eGFP)
avg <- AverageExpression(data_eGFP, features = hvg, group.by = group_by, assays = "RNA", slot = "data")$RNA
avg <- avg[apply(avg, 1, var) > 0, , drop = FALSE]

# 2) seleccionar top-N por varianza (opcional)
if (!is.null(top_var_rows) && top_var_rows < nrow(avg)) {
  v <- apply(avg, 1, var)
  avg <- avg[order(v, decreasing = TRUE)[seq_len(top_var_rows)], , drop = FALSE]
}

# 3) z-score por fila y capping
z <- t(scale(t(as.matrix(avg))))
z[z < cap_limits[1]] <- cap_limits[1]
z[z > cap_limits[2]] <- cap_limits[2]

# 4) función de color
make_col_fun <- function(low, high, palette=c("viridis","RdBu","PuOr")[2]) {
  if (palette=="viridis") cols <- viridisLite::viridis(11) else
  if (palette=="RdBu")   cols <- rev(brewer.pal(11,"RdBu")) else
                           cols <- colorRampPalette(brewer.pal(9,"PuOr"))(11)
  circlize::colorRamp2(seq(low, high, length.out=length(cols)), cols)
}
col_fun <- make_col_fun(cap_limits[1], cap_limits[2], palette_cont)

# 5) distancia tipo 1 - correlación (mejor para perfiles de expresión)
dist_cols <- as.dist(1 - cor(z, method="pearson"))
dist_rows <- as.dist(1 - cor(t(z), method="pearson"))

ht <- Heatmap(
  z, name="Z-score", col=col_fun,
  clustering_distance_columns = dist_cols,
  clustering_distance_rows    = dist_rows,
  clustering_method_columns   = "average",
  clustering_method_rows      = "average",
  show_row_names = FALSE, show_column_names = TRUE,
  column_title = "Average expression (z-score per gene)",
  border=TRUE, rect_gp = gpar(col="#EEEEEE", lwd=0.4),
  heatmap_legend_param = list(at=seq(cap_limits[1],cap_limits[2],by=1),
                              title_gp=gpar(fontface="bold"), labels_gp=gpar(fontsize=8)),
  use_raster=TRUE, raster_quality=3
)

# 6) exportación vectorial
dims_w <- 3.5 + 0.22 * ncol(z); dims_h <- 4.0 + 0.18 * nrow(z)
dims_w <- max(6, min(dims_w, 16)); dims_h <- max(6, min(dims_h, 20))

grDevices::cairo_pdf(outfile, width=5, height=10, family="Arial")
draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
message("Guardado: ", outfile)
