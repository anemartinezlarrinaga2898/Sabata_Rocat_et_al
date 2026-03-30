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
library(optparse)
library(pheatmap)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Correlation",sep="/")
dir.create(path.guardar,recursive=TRUE)

ident<-"AnnotLayer"

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))
################################################################################

# --- Carga y preparación ---

# | WILCOXON: 
path_wilcoxon<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Wilcoxon/eGFP_WT_MUT",sep="/")
file_name_wilcoxon<-paste("DEG_Mut_WT_eGFP_Pos_",ident,".rds",sep="")

markers_wilcoxon<-readRDS(file.path(path_wilcoxon,file_name_wilcoxon))
markers_wilcoxon<-do.call(rbind,markers_wilcoxon)

# | DESEQ2: 
path_deseq<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Pseudobulk/eGFP_WT_MUT",sep="/")
list_cell_types<-list.files(path_deseq)
file_name_deseq2<-"Res_Repro.rds"

Lista_Deseq_Values<-list()
for(i in seq_along(list_cell_types)){
    ct<-list_cell_types[i]
    path_file<-paste(path_deseq,ct,sep="/")
    markers_deseq2<-readRDS(file.path(path_file,file_name_deseq2))
    markers_deseq2$cluster<-ct
    Lista_Deseq_Values[[i]]<-markers_deseq2
}

markers_deseq2<-as.data.frame(do.call(rbind,Lista_Deseq_Values))

# | Joined datasets 

wilc <- markers_wilcoxon %>% dplyr::select(cluster, gene, avg_log2FC)
deseq <- markers_deseq2 %>% dplyr::select(cluster, gene, log2FoldChange)

df_all <- inner_join(wilc, deseq,by = c("cluster", "gene"))

# | Correlation per cluster 

cor_table <- df_all %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    spearman_cor = cor(avg_log2FC,
                       log2FoldChange,
                       method = "spearman"),
    pearson_cor  = cor(avg_log2FC,
                       log2FoldChange,
                       method = "pearson"),
    n_genes = n())

# ..................................................
# Heatmap: Correlation Spearson 
# ..................................................

mat_cor <- cor_table %>%
  select(cluster, spearman_cor) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

df_plot <- as.data.frame(mat_cor) %>%
  rownames_to_column("Cluster")

ggplot(df_plot, aes(x = "Correlation", y = Cluster, fill = spearman_cor)) +
   geom_tile(color = "white", size = 0.8) +
   geom_text(aes(label = round(spearman_cor, 3)), size = 4)+
  scale_fill_gradient(
    low = "#E3F2FD",
    high = "#0D47A1",
    limits = c(0.85, 0.91)
  ) +
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
  labs(x = "", y = "")
ggsave(file.path(path.guardar,"Correlation_Spearman.png"),width=3,heigh=5)
ggsave(file.path(path.guardar,"Correlation_Spearman.pdf"),width=3,heigh=5)
# ..................................................
# Heatmap: Correlation Spearson 
# ..................................................

mat_cor <- cor_table %>%
  select(cluster, pearson_cor) %>%
  tibble::column_to_rownames("cluster") %>%
  as.matrix()

df_plot <- as.data.frame(mat_cor) %>%
  rownames_to_column("Cluster")

ggplot(df_plot, aes(x = "Correlation", y = Cluster, fill = pearson_cor)) +
   geom_tile(color = "white", size = 0.8) +
   geom_text(aes(label = round(pearson_cor, 3)), size = 4)+
   scale_fill_gradient(
    low = "#E0F2F1",
    high = "#004D40",
    limits = c(0.6, 0.91)) +
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
  labs(x = "", y = "")
ggsave(file.path(path.guardar,"Correlation_Pearson.png"),width=3,heigh=5)
ggsave(file.path(path.guardar,"Correlation_Pearson.pdf"),width=3,heigh=5)

# ..................................................
# Scatter plots per cell types 
# ..................................................

ggplot(df_all, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~cluster) +
  theme_bw()
ggsave(file.path(path.guardar,"ScatterPlots.png"),width=5,heigh=5)
ggsave(file.path(path.guardar,"ScatterPlots.pdf"),width=5,heigh=5)

# ..................................................
# Porcentaje de concordancia de dirección
# ..................................................

direction_summary <- df_all %>%
  mutate(
    dir_w = sign(avg_log2FC),
    dir_d = sign(log2FoldChange),
    concordant = dir_w == dir_d) %>%
  group_by(cluster) %>%
  summarise(
    percent_concordant = mean(concordant) * 100)

ggplot(direction_summary,
       aes(x = cluster, y = percent_concordant)) +
  geom_col(fill = "#1976D2", width = 0.6) +
  geom_text(aes(label = paste0(round(percent_concordant,1), "%")),
            vjust = -0.5,
            fontface = "bold") +
  ylim(0,100) +
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
  labs(title = "",
       y = "% genes with same direction",
       x = "")
ggsave(file.path(path.guardar,"BarPlot_Direction.png"),width=5,heigh=5)
ggsave(file.path(path.guardar,"BarPlot_Direction.pdf"),width=5,heigh=5)

# ..................................................
# | Estimate sig genes: 
# ..................................................

sig_w <- markers_wilcoxon %>%
  filter(p_val_adj < 0.05)

sig_d <- markers_deseq2 %>%
  filter(padj < 0.05)

overlap_summary <- sig_w %>%
  inner_join(sig_d, by = c("cluster", "gene")) %>%
  group_by(cluster) %>%
  summarise(n_overlap = n())

count_w <- sig_w %>%
  group_by(cluster) %>%
  summarise(n_wilcoxon = n())

count_d <- sig_d %>%
  group_by(cluster) %>%
  summarise(n_deseq = n())

overlap_table <- count_w %>%
  left_join(count_d, by="cluster") %>%
  left_join(overlap_summary, by="cluster") %>%
  mutate(
    percent_overlap_w = n_overlap / n_wilcoxon * 100,
    percent_overlap_d = n_overlap / n_deseq * 100)


overlap_long <- overlap_table %>%
  select(cluster, n_wilcoxon, n_overlap) %>%
  pivot_longer(
    cols = c(n_wilcoxon, n_overlap),
    names_to = "Method",
    values_to = "Count"
  ) %>%
  mutate(
    Method = recode(Method,
                    n_wilcoxon = "Wilcoxon (significant)",
                    n_overlap  = "Overlap with DESeq2"))

ggplot(overlap_long,
       aes(x = cluster, y = Count, fill = Method)) +
  geom_col(position = "identity", width = 0.6) +
  scale_fill_manual(
    values = c("Wilcoxon (significant)" = "#90CAF9",
               "Overlap with DESeq2"   = "#1565C0")) +
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
  labs(
    title = "Significant gene overlap",
    y = "Number of genes",
    x = "",
    fill = "")
ggsave(file.path(path.guardar,"BarPlot_SigGenesOverlap.png"),width=6,heigh=4)
ggsave(file.path(path.guardar,"BarPlot_SigGenesOverlap.pdf"),width=6,heigh=4)

# ..................................................
# | Genes Clave
# ..................................................

genes_check<-c("Gja4", "Unc5b", "Cxcr4", "Esm1", "Hey1", "Sat1","Tll1","Flrt2","Lama1","Nr2f2","Tgfbi","Mafb","Adgrg6")

# --- Wilcoxon ---
wilc_long <- markers_wilcoxon %>%
  transmute(
    cluster = cluster,
    gene    = gene,
    method  = "Wilcoxon",
    log2FC  = avg_log2FC,
    padj    = p_val_adj)

# --- DESeq2 ---
deseq_long <- markers_deseq2 %>%
  transmute(
    cluster = cluster,
    gene    = gene,            # si en tu objeto es SYMBOL: gene = SYMBOL
    method  = "DESeq2",
    log2FC  = log2FoldChange,
    padj    = padj)

df_methods <- bind_rows(wilc_long, deseq_long) %>%
  filter(!is.na(gene), !is.na(cluster), !is.na(log2FC))
df_methods_filter<-df_methods[which(df_methods$gene %in% genes_check),]

df_dot <- df_methods_filter %>%
  mutate(
    neglog10_padj = -log10(padj + 1e-300),
    signif = ifelse(padj < 0.05, "Significant", "Not Significant"))

ggplot(df_dot, aes(x = method, y = gene, color = log2FC,size=neglog10_padj)) +
  geom_point() +
  facet_grid(cols=vars(cluster),scales="free") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(color = "log2FC")
ggsave(file.path(path.guardar,"DotPlot_Specific_Genes.png"),width=8.5,heigh=4)
ggsave(file.path(path.guardar,"DotPlot_Specific_Genes.pdf"),width=8.5,heigh=4)
