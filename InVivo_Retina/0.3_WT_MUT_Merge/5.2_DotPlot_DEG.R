################################################################################
# SCRIPT: Visualization of selected WT versus MUT differential-expression results
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Builds heatmap and dot-plot summaries from previously calculated differential-expression results and a curated gene/signature selection.
#
# INPUT:
#   results/0.6_Join_MUT_WT/DownStream/DEG/Annotations/DEG_Mut_WT_eGFP_Pos_Layer_1.rds
#   results/0.6_Join_MUT_WT/DownStream/DEG/Annotations/DEG_Mut_WT_eGFP_Pos_Layer_1.xlsx (sheet: Opcion1)
#
# OUTPUT:
#   Selected DEG heatmap and signature DotPlots
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

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/DownStream/DEG/Annotations/Opcion2",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# Load project utility functions 
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

args <- commandArgs(trailingOnly = TRUE)
################################################################################

# Data loading and preparation
print("Load data")
# Path obj: 
path.guardar_original <- results_root
path.obj<-paste(path.guardar_original,"0.6_Join_MUT_WT/DownStream/DEG/Annotations",sep="/")
ident<-"Layer_1"
file_name <- paste("DEG_Mut_WT_eGFP_Pos_", ident, ".rds", sep = "")
lista_dfs<-readRDS(file.path(path.obj, file_name))
markers_logFC<-do.call(rbind,lista_dfs)

markers<-readxl::read_excel(file.path(results_root, "0.6_Join_MUT_WT", "DownStream", "DEG", "Annotations", "DEG_Mut_WT_eGFP_Pos_Layer_1.xlsx"),sheet="Opcion1")

# Select the genes of interest
markers_logFC_filter<-markers_logFC[which(markers_logFC$gene %in% markers$Gene),]

# --------------------------------------------------
# Heatmap LogFC 
# --------------------------------------------------

mat <- markers_logFC_filter %>%
  select(gene, cluster, avg_log2FC) %>%
  distinct(gene, cluster, .keep_all = TRUE) %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# recorte opcional (winsorize) para que el color sea más informativo
mat2 <- pmax(pmin(mat, 2), -2)

pdf(paste(path.guardar, "Prueba.pdf", sep = "/"), width = 5, height = 10) 
# Generate the heatmap with the specified parameters
Heatmap(mat, name = "avg_log2FC")
dev.off()

# --------------------------------------------------
# DotPlot of the logFC
# --------------------------------------------------

df_dot <- markers_logFC_filter %>%
  mutate(
    neglog10_padj = -log10(p_val_adj + 1e-300)  # evita Inf si hay 0
  )

ggplot(df_dot, aes(x = cluster, y = gene)) +
  geom_point(aes(size = neglog10_padj, color = avg_log2FC), alpha = 0.9) +
  scale_size_continuous(name = "-log10(p_adj)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold")
  ) +
  labs(x = "Cluster", y = "Gene", color = "avg_log2FC")+
   scale_color_gradient2(
    low = "#2166AC",     # azul
    mid = "white",
    high = "#B2182B",    # rojo
    midpoint = 0,
    name = "avg_log2FC"
  ) 
ggsave(file=paste(path.guardar,"DotPlot_PI3K_Signatures.pdf",sep="/"),width=5,height=10)

# ------------------------------------------------------------------------------------------------
# DotPlot: Chaning the size so the signicant is bigger and the not significant is less 
# ------------------------------------------------------------------------------------------------
df_dot <- markers_logFC_filter %>%
  mutate(
    neglog10_padj = -log10(p_val_adj + 1e-300),
    signif = ifelse(p_val_adj < 0.05, "Significant", "Not Significant")
  )

ggplot(df_dot, aes(x = cluster, y = gene)) +
  geom_point(aes(size = signif, color = avg_log2FC), alpha = 0.9) +
  scale_size_manual(
    values = c("Not Significant" = 1, "Significant" = 5),
    name = "Significance"
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "avg_log2FC"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold")
  ) +
  labs(x = "Cluster", y = "Gene")

ggsave(file=paste(path.guardar,"DotPlot_PI3K_Signatures_Size_PValue_BigSmall.pdf",sep="/"),width=5,height=12)

# ------------------------------------------------------------------------------------------------
# Dotplot one per signature
# ------------------------------------------------------------------------------------------------

markers$Features<-stringr::str_replace_all(markers$Features,"/","_")
markers$Features<-stringr::str_replace_all(markers$Features,"-","_")
markers$Features<-stringr::str_replace_all(markers$Features," ","_")
signatures<-unique(markers$Features)

for(i in seq_along(signatures)){
  sig<-signatures[i]
  print(sig)
  markers_sig<-markers[which(markers$Features==sig),]
  markers_logFC_filter<-markers_logFC[which(markers_logFC$gene %in% markers_sig$Gene),]

  df_dot <- markers_logFC_filter %>%
  mutate(
    neglog10_padj = -log10(p_val_adj + 1e-300),
    signif = ifelse(p_val_adj < 0.05, "Significant", "Not Significant"))
  gene_levels <- unique(markers_sig$Gene)
  df_dot$gene <- factor(df_dot$gene, levels = gene_levels)

  ggplot(df_dot, aes(x = cluster, y = gene)) +
    geom_point(aes(size = signif, color = avg_log2FC), alpha = 0.9) +
    scale_size_manual(values = c("Not Significant" = 1, "Significant" = 5),name = "Significance") +
    scale_color_gradient2(low = "#2166AC",mid = "white",high = "#B2182B",midpoint = 0,name = "avg_log2FC") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold")) +
    labs(x = "Cluster", y = "Gene")
  file_name<-paste("DotPlot_Signatures_logFC_PValue_Size_",sig,".pdf",sep="")
  ggsave(file=paste(path.guardar,file_name,sep="/"),width=5,height=5)
}

