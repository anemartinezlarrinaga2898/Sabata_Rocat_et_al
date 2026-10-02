################################################################################
# SCRIPT: eGFP distribution after MUT endothelial cluster removal
# AUTHOR: Ane Martinez Larrinaga
# DATE: 19-12-2023
#
# DESCRIPTION:
# Defines eGFP-positive cells using the original expression threshold, visualizes eGFP distribution after cluster removal, and retains the exploratory cell-cycle and pre-arterial signature analyses present in the original script.
#
# INPUT:
#   results/0.5_MUT_Analysis/EC_Subset/RemoveClusters/Remove_Res_05_Clus3679/Harmony/Harmony.rds
#
# OUTPUT:
#   eGFP distribution and exploratory analysis plots.
#
# NOTE:
#   Analytical thresholds, clustering resolutions and biological selections
#   are preserved from the original analysis unless explicitly documented in
#   REVIEW_NOTES.md.
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

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(readxl)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.5_MUT_Analysis/EC_Subset/RemoveClusters/Remove_Res_05_Clus3679/Harmony/egfp_Levels",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load project utility functions
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))
source(file.path(utils_dir, "0.19_Util_CellType_Classifier.R"))
# ccAFv2 color palette
cols_ccaf <- c('G1' = '#f37f73', 'G2/M' = '#3db270', 'Late G1' = '#1fb1a9',
               'M/Early G1' = '#6d90ca', 'Neural G0' = '#d9a428', 'S' = '#8571b2', 
               'S/G2' = '#db7092', 'G0/G1' = '#FF6600', 'Unknown' = '#cccccc')
args <- commandArgs(trailingOnly = TRUE)
################################################################################

path.obj<-paste(results_root,"0.5_MUT_Analysis/EC_Subset/RemoveClusters/Remove_Res_05_Clus3679/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"
if (length(args) < 1) {
  stop(
    "Missing clustering identity. Run as: Rscript 1.5_eGFP_Dist.R <metadata_column>"
  )
}

ident <- args[1]

if (!ident %in% colnames(data@meta.data)) {
  stop(paste("Metadata column not found:", ident))
}


data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == "TRUE" ~ "eGFP +",eGFP_pos == "FALSE" ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# Visualize eGFP distribution

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster = FALSE)&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_Total_eGFP_Levels_",ident,".png",sep="_")),width = 6,height = 6)

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster = FALSE,split.by="Phenotype")&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_ByPheno_eGFP_Levels_",ident,".png",sep="_")),width = 12,height = 6)

# Quantify eGFP distribution

# Extract metadata
df_meta <- data@meta.data %>% 
  select(Cluster = all_of(ident), Phenotype, eGFP_Levels) %>% 
  filter(!is.na(eGFP_Levels)) # Remove missing classifications

# Calculate eGFP proportions by cluster and phenotype
df_porcentajes <- df_meta %>%
  group_by(Cluster, Phenotype, eGFP_Levels) %>%
  tally() %>%
  group_by(Cluster, Phenotype) %>%
  mutate(Porcentaje = (n / sum(n)) * 100)

# Plot eGFP proportions
ggplot(df_porcentajes, aes(x = Cluster, y = Porcentaje, fill = eGFP_Levels)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Phenotype, ncol = 2) + # Separate phenotypes
  scale_fill_manual(values = c("eGFP +" = "blue", "eGFP -" = "grey80")) + # eGFP-level colors
  labs(y = "Porcentaje de células (%)", x = "Clusters", fill = "eGFP Levels") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "plain"),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"), # Estilo del título de la faceta (WT/MUT)
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

# Save bar plot
file_name_bar <- paste("BarPlot_eGFP_Levels_by_Cluster_Pheno_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

# Marker DotPlot

# Marker genes used for cluster visualization
genes_dotplot <- c(
  # Cluster 0: Pool Proliferativo (Mitosis residual alta)
  "Top2a", "Mki67", 
  # Cluster 1: Barrera Hematorretiniana / Maduración metabólica
  "Apod", "Slc38a5", "Trf",
  # Cluster 2: Tip Cells de Vanguardia (Hiper-angiogénicas)
  "Kcne3", "Angpt2", "Esm1",
  # Cluster 3: Transición Tip-to-Artery (Zonación temprana)
  "Bmx", "Gja5", "Dkk2",
  # Cluster 4: Transición Venosa / Remodelación de tallo
  "Nr2f2", "Bgn",
  # Cluster 5: Estado Inmuno-Endotelial (Respuesta a IFN)
  "Six3", "Ifi44", "Cfh"
)

# Retain only genes present in the RNA assay
genes_dotplot <- intersect(genes_dotplot, rownames(data))

# Generate DotPlot
DotPlot(data, features = genes_dotplot, assay = "RNA",group.by=ident) &
  coord_flip() & # Voltea el gráfico para que los genes queden en el eje Y (más legible)
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "italic"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) # Paleta clásica azul/rojo muy limpia

# Save DotPlot
file_name_dot <- paste("DotPlot_Identity_Markers_", ident, ".pdf", sep = "")
ggsave(filename = file.path(path.guardar, file_name_dot), width = 8, height = 7, dpi = 300)

# Phenotype distribution

matchSCore2::summary_barplot(class.fac = data$Phenotype,obs.fac =data@active.ident)
file_name_bar <- paste("BarPlot_Phenotype_Distribution_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

DimPlot(data,reduction = "umap",group.by = ident,label = TRUE, label.size = 5,pt.size = 1.5,raster = FALSE)&NoAxes()&NoLegend()
file_name_bar <- paste("Dimplot_Clusters_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

DimPlot(data,reduction = "umap",group.by = ident,label = TRUE, label.size = 5,pt.size = 1.5,raster = FALSE,split.by="Phenotype")&NoAxes()&NoLegend()
file_name_bar <- paste("Dimplot_Clusters_ByPheno", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)


# WT subset used in the original exploratory analysis

data<-SetIdent(data,value="Phenotype")
data_wt<-subset(data,ident="WT")

# Estimate cell-cycle states

data_wt <- PredictCellCycle(data_wt,
                 threshold=0.5,
                 include_g0 = TRUE,
                 do_sctransform=TRUE,
                 assay='RNA',
                 species='mouse',
                 gene_id='symbol',
                 spatial = FALSE)

# Cell-cycle state UMAP
DimPlot(data_wt, group.by = "ccAFv2",cols=cols_ccaf,pt.size = 1.5)&NoAxes()
file_name_bar <- paste("Dimplot_Stages_ccaFv2_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

matchSCore2::summary_barplot(class.fac = data_wt$ccAFv2,obs.fac =data_wt@active.ident) +scale_fill_manual(values = cols_ccaf) 
file_name_bar <- paste("BarPlot_ccaFv2_Distribution_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

# Cell-cycle composition by cluster and phenotype
df_prop <- data_wt@meta.data %>%
  dplyr::group_by(.data[[ident]], Phenotype, ccAFv2) %>%
  tally() %>%
  group_by(.data[[ident]], Phenotype) %>%
  mutate(pct = n / sum(n) * 100)

ggplot(df_prop, aes(x = Phenotype, y = pct, fill = ccAFv2)) +
  geom_bar(stat = "identity") +
 facet_wrap(vars(.data[[ident]])) +
  scale_fill_manual(values = cols_ccaf) +
  labs(y = "Porcentaje de células", title = "Composición del ciclo por Layer")
file_name_bar <- paste("BarPlot_ccaFv2_by_Pheno_In_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

# Visualize selected cell-cycle scores
FeaturePlot(data_wt, features = "Late.G1", order = TRUE,label = TRUE)
file_name_bar <- paste("Fp_Late_G1_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

FeaturePlot(data_wt, features = "M.Early.G1", order = TRUE,label = TRUE)
file_name_bar <- paste("Fp_Early_G1_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

VlnPlot(data_wt,features = "Late.G1",sort = TRUE)
file_name_bar <- paste("Vln_Late_G1_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

VlnPlot(data_wt,features = "M.Early.G1",sort = TRUE)
file_name_bar <- paste("Vln_Early_G1_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

# Estimate pre-arterial signature

genes_pre_arteriales <- c("Unc5b", "Dll4", "Cxcr4") # Pon tus genes aquí
data_wt <- UCell::AddModuleScore_UCell(data_wt, features = list(PreArterial = genes_pre_arteriales), name = NULL)

# Visualize selected cell-cycle scores
FeaturePlot(data_wt, features = "PreArterial", order = TRUE,label = TRUE)
file_name_bar <- paste("Fp_PreArterial_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

VlnPlot(data_wt,features = "PreArterial",sort = TRUE)
file_name_bar <- paste("Vln_PreArterial_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

# Correlate pre-arterial signature with cell-cycle scores

df_corr <- FetchData(data_wt, vars = c(ident, "PreArterial", "M.Early.G1","Late.G1","G2.M","S.G2","S","G1", "Phenotype")) %>%
  filter(!is.na(.data[[ident]]))

cor_results <- df_corr %>%
  dplyr::group_by(.data[[ident]]) %>%
  dplyr::summarize(
    Corr_M.Early.G1 = cor(PreArterial, M.Early.G1, method = "spearman"),
    Corr_Late.G1  = cor(PreArterial, Late.G1, method = "spearman"),
    Corr_G2.M  = cor(PreArterial, G2.M, method = "spearman"),
    Corr_S.G2  = cor(PreArterial, S.G2, method = "spearman"),
    Corr_S  = cor(PreArterial, S, method = "spearman"),
    Corr_G1  = cor(PreArterial, G1, method = "spearman")
  ) %>%
  pivot_longer(cols = starts_with("Corr"), names_to = "Phase", values_to = "Correlation")

ggplot(cor_results, aes(x = Phase, y = .data[[ident]], fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0) +
  labs(title = "Correlación: Identidad Pre-Arterial vs Ciclo Celular",
       subtitle = "Coeficiente de Spearman por Cluster",
       fill = "Rho")+
  theme_classic() +
  theme(
    # angle = 90 rota el texto
    # vjust controla la posición vertical
    # hjust = 1 alinea el texto al eje para que no flote
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
file_name_bar <- paste("Correlation_ccaFv2_early_late_prearterialsig_WT_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)
