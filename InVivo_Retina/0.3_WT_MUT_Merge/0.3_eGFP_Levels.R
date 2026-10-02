################################################################################
# SCRIPT: eGFP distribution in the combined WT + MUT dataset
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Defines eGFP-positive cells using the original expression threshold, visualizes eGFP distribution by phenotype and cluster, and generates marker and composition plots.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
#
# OUTPUT:
#   eGFP distribution, cluster-composition and marker plots
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
library(readxl)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Harmony/egfp_Levels",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# Load project utility functions 
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))
source(file.path(utils_dir, "0.19_Util_CellType_Classifier.R"))
# Standard ccAFv2 color palette
cols_ccaf <- c('G1' = '#f37f73', 'G2/M' = '#3db270', 'Late G1' = '#1fb1a9',
               'M/Early G1' = '#6d90ca', 'Neural G0' = '#d9a428', 'S' = '#8571b2', 
               'S/G2' = '#db7092', 'G0/G1' = '#FF6600', 'Unknown' = '#cccccc')
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !nzchar(args[1])) {
  stop("Please provide the metadata identity column as the first command-line argument.")
}
################################################################################

path.obj<-paste(results_root,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"
ident<-args[1]

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# | Dimplot of EGFP Distribution 

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_Total_eGFP_Levels_",ident,".png",sep="_")),width = 6,height = 6)

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster=FALSE,split.by="Phenotype")&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_ByPheno_eGFP_Levels_",ident,".png",sep="_")),width = 12,height = 6)


# | Dimplot of EGFP Distribution 

data <- SetIdent(data, value = "eGFP_pos")
data_egfp <- subset(data, idents = "TRUE")

# 1. Extraer metadatos a un dataframe
df_meta <- data_egfp@meta.data %>% 
  select(Cluster = all_of(ident), Phenotype)

# 2. Calcular porcentajes por Cluster, Fenotipo y Nivel de eGFP
df_porcentajes <- df_meta %>%
  group_by(Cluster, Phenotype) %>%
  tally() %>%
  group_by(Cluster, Phenotype) %>%
  mutate(Porcentaje = (n / sum(n)) * 100)

# 3. Graficar con ggplot2
ggplot(df_porcentajes, aes(x = Phenotype, y = Porcentaje, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(y = "", x = "", fill = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "plain"),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"), # Estilo del título de la faceta (WT/MUT)
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

# 4. Guardar el Barplot
file_name_bar <- paste("BarPlot_eGFP_Levels_by_Cluster_Pheno_OnlyEGFP_Pos_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)


# | BarPlot Distirbution 

# 1. Extraer metadatos a un dataframe
df_meta <- data@meta.data %>% 
  select(Cluster = all_of(ident), Phenotype, eGFP_Levels) %>% 
  filter(!is.na(eGFP_Levels)) # Por si acaso hay NAs

# 2. Calcular porcentajes por Cluster, Fenotipo y Nivel de eGFP
df_porcentajes <- df_meta %>%
  group_by(Cluster, Phenotype, eGFP_Levels) %>%
  tally() %>%
  group_by(Cluster, Phenotype) %>%
  mutate(Porcentaje = (n / sum(n)) * 100)

# 3. Graficar con ggplot2
ggplot(df_porcentajes, aes(x = Cluster, y = Porcentaje, fill = eGFP_Levels)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Phenotype, ncol = 2) + # Separa WT y MUT
  scale_fill_manual(values = c("eGFP +" = "blue", "eGFP -" = "grey80")) + # Ajusta según tus nombres exactos
  labs(y = "Porcentaje de células (%)", x = "Clusters", fill = "eGFP Levels") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "plain"),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 13, face = "bold"), # Estilo del título de la faceta (WT/MUT)
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

# 4. Guardar el Barplot
file_name_bar <- paste("BarPlot_eGFP_Levels_by_Cluster_Pheno_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

# | Dot Plots markers 

# Definimos una lista limpia con los mejores marcadores de identidad por cluster
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

# Asegurar que los genes seleccionados existen en la matriz de RNA
genes_dotplot <- intersect(genes_dotplot, rownames(data))

# Generar el DotPlot estructurado
DotPlot(data, features = genes_dotplot, assay = "RNA",group.by=ident) &
  coord_flip() & # Voltea el gráfico para que los genes queden en el eje Y (más legible)
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "italic"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) # Paleta clásica azul/rojo muy limpia

# Guardar el DotPlot
file_name_dot <- paste("DotPlot_Identity_Markers_", ident, ".pdf", sep = "")
ggsave(filename = file.path(path.guardar, file_name_dot), width = 8, height = 7, dpi = 300)

# | BarPlot Phenotype

matchSCore2::summary_barplot(class.fac = data$Phenotype,obs.fac =data@active.ident)
file_name_bar <- paste("BarPlot_Phenotype_Distribution_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)

DimPlot(data,reduction = "umap",group.by = ident,label = TRUE, label.size = 5,pt.size = 1.5,raster=FALSE)&NoAxes()&NoLegend()
file_name_bar <- paste("Dimplot_Clusters_", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 5, height = 5, dpi = 300)

DimPlot(data,reduction = "umap",group.by = ident,label = TRUE, label.size = 5,pt.size = 1.5,raster=FALSE,split.by="Phenotype")&NoAxes()&NoLegend()
file_name_bar <- paste("Dimplot_Clusters_ByPheno", ident, ".png", sep = "")
ggsave(filename = file.path(path.guardar, file_name_bar), width = 10, height = 5, dpi = 300)


