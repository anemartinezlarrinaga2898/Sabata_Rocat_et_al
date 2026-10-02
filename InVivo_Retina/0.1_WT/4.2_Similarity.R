################################################################################
# PUBLICATION CODE
# FILE: 5.2_Similarity.R
# PURPOSE: Assess similarity among final annotated endothelial populations.
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
library(SingleCellExperiment)
library(ggpubr)
library(ggplot2)
library(decoupleR)
library(OmnipathR)
library(viper)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(dendextend)
library(pheatmap)
library(ggdendro)
library(pvclust)
library(plyr)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))


# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)
################################################################################

path.guardar_original <- PATHS$results
path.obj<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
data <- readRDS(paste(path.obj,"Data_WT_Anotado.rds",sep="/"))
data<-SetIdent(data,value="AnnotLayer")
print("data loaded")

args = commandArgs(trailingOnly=TRUE)
ident<-"AnnotLayer"

data<-SetIdent(data,value=ident)

# 2) Compute average expression matrix (genes × tissues)
#    Using only highly variable genes (HVGs) from your Seurat object
hvg_genes <- VariableFeatures(data)  
expr_mat  <- AverageExpression(
  data,
  features = hvg_genes,
  group.by  = ident,
  assays    = "RNA",
  slot      = "data"
)$RNA
scaled_matrix <- t(scale(t(expr_mat)))

# Remove any genes with zero variance
expr_mat <- expr_mat[apply(expr_mat, 1, var) > 0, ]
expr_mat_dense <- as.matrix(expr_mat)

# 3) Hierarchical clustering with multiscale bootstrap (pvclust)
set.seed(42)  # for reproducibility
pv <- pvclust(
  expr_mat_dense,
  method.hclust = "ward.D2",  # clustering linkage
  method.dist   = "euclidean",# distance metric
  nboot         = 1000        # number of bootstrap replicates
)
file_name<-paste(ident,"pv.rds",sep="_")
saveRDS(pv,paste(path.guardar,file_name,sep="/"))

# Plot dendrogram showing AU (red) and BP (green) p-values
file_name<-paste(ident,"Dend.pdf",sep="_")
pdf(paste(path.guardar,file_name,sep="/"),width=5,height=5)
plot(pv)
# Highlight clusters with AU ≥ 0.95
pvrect(pv, alpha = 0.95)
dev.off()

# 4) Correlation heatmap of tissue signatures
# 4.1) Calculate Pearson correlation between tissues
cor_mat <- cor(expr_mat_dense, method = "pearson")

# 4.2) Option A: pheatmap
file_name<-paste(ident,"Heatmap_Correlation.pdf",sep="_")
palette_pheat <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
pdf(paste(path.guardar,file_name,sep="/"),width=5,height=5)
pheatmap(
  cor_mat,
  color                   = palette_pheat,
  clustering_method       = "ward.D2",
  clustering_distance_rows   = "euclidean",
  clustering_distance_cols   = "euclidean",
  treeheight_row          = 40,
  treeheight_col          = 40,
  border_color            = NA,
  main                    = "Pearson correlation of endothelial populations"
)
dev.off()

########################################################################################################################################

# --- PARÁMETROS QUE PUEDES EDITAR ---
file_name<-paste(ident,"HC_AverageExpression.pdf",sep="_")
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
file_name<-paste(ident,"Matrix_Heatmap.rds",sep="_")
saveRDS(z,file.path(path.guardar,file_name))
file_name<-paste(ident,"Matrix_Heatmap_HC.xlsx",sep="_")
openxlsx::write.xlsx(z,paste(path.guardar,file_name,sep="/"))


# 6) exportación vectorial
dims_w <- 3.5 + 0.22 * ncol(z); dims_h <- 4.0 + 0.18 * nrow(z)
dims_w <- max(6, min(dims_w, 16)); dims_h <- max(6, min(dims_h, 20))

grDevices::cairo_pdf(outfile, width=5, height=10, family="Arial")
draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()
message("Guardado: ", outfile)
