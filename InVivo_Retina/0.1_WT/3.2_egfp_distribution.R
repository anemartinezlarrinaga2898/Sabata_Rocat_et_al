################################################################################
# PUBLICATION CODE
# FILE: 3.2_egfp_distribution.R
# PURPOSE: Assess eGFP-positive cell distribution across clusters and phenotypes.
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
library(readxl)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239/egfp_Levels",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")
source_required("0.19_Util_CellType_Classifier.R")
args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-args[1]

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# | Dimplot of EGFP Distribution 

DimPlot(data,reduction = "umap",group.by = ident,label = FALSE, label.size = 5,pt.size = 1.5,raster=FALSE)&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_Total_Clusters_",ident,".png",sep="_")),width = 6,height = 6)

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_Total_eGFP_Levels_",ident,".png",sep="_")),width = 6,height = 6)

DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster=FALSE,split.by="Phenotype")&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_ByPheno_eGFP_Levels_",ident,".png",sep="_")),width = 12,height = 6)

# | BarPlot Distirbution 

data<-SetIdent(data,value=ident)

matchSCore2::summary_barplot(class.fac = data$Phenotype,obs.fac =data@active.ident)
ggsave(file.path(path.guardar,paste("BarPlot_Pheno_",ident,".png",sep="_")),width=5,height=5)

matchSCore2::summary_barplot(class.fac = data@active.ident,obs.fac =data$Phenotype)
ggsave(file.path(path.guardar,paste("BarPlot_",ident,"_Pheno",".png",sep="_")),width=5,height=5)

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

