################################################################################
# SCRIPT: Initial annotation of the combined WT + MUT endothelial dataset
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Assigns endothelial population labels to the joint WT + MUT dataset using the original manual cluster mapping, defines eGFP-positive cells, and generates annotation and marker visualizations.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Harmony/Harmony.rds
#   data/WT_MUT_markers.xlsx (sheet: Hoja2)
#
# OUTPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Annotations/Data_EC_Annotated.rds and annotation figures
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
library(plyr)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Annotations",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
################################################################################

path.obj<-paste(results_root,"0.6_Join_MUT_WT/EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-"Harmony_Log_res.0.9"
data@meta.data$Layer_1<- revalue(data@meta.data[[ident]], c("0" = "VenousCapillary EC",
                                                               "1" = "Proliferative 1",
                                                               "2" = "Proliferative 2",
                                                               "3" = "VenousCapillary EC",
                                                               "4" = "Angiogenic EC",
                                                               "5" = "ArteryCapillary EC",
                                                               "6" = "ArteryCapillary EC",
                                                               "7" = "Proliferative 1",
                                                               "8" = "Angiogenic EC",
                                                               "9" = "VenousCapillary EC",
                                                               "10" = "Arterial EC"))

idx_tip<-which(data$Harmony_Log_res.0.5=="8")
data@meta.data$Layer_1<-as.character(data@meta.data$Layer_1)
data@meta.data$Layer_1[idx_tip]<-"Tip EC"

saveRDS(data,file.path(path.guardar,"Data_EC_Annotated.rds"))

ident<-"Layer_1"
data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"))

# Define colors: 
# ----------------------------------------
# UMAP All populations
cluster_order<-c("Tip EC","Angiogenic EC","Proliferative 1","Proliferative 2","VenousCapillary EC","ArteryCapillary EC","Arterial EC")
colors_use<-c("#b37087","#f4cedb","#b5dcfb","#90caf9","#C5E1A5","#FFE082","#FF8A65")
names(colors_use) <- cluster_order
# ----------------------------------------

data$Layer_1<-factor(data$Layer_1,levels=cluster_order)

# UMAP Total
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = TRUE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_Label.pdf",sep="/"),width=7,height=7)
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = FALSE, cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabel.pdf",sep="/"),width=7,height=7)
DimPlot(data,reduction = "umap",group.by = "Layer_1",label = FALSE, cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabelNoLegend.pdf",sep="/"),width=7,height=7)

# ----------------------------------------
# Violin Plot of EC Markers 
# ----------------------------------------

# 2. Definir los genes de interés (formato correcto: Primera mayúscula para ratón)
genes_endoteliales <- c("Pecam1", "Cdh5")

# 3. Generar el VlnPlot
vln_plot <- VlnPlot(
  data, 
  features = genes_endoteliales, 
  group.by = "Layer_1", 
  cols = colors_use,       # Aplicamos tu paleta homogenizada
  pt.size = 0.1,                # Tamaño de los puntos (células) sobre el violín
  combine = TRUE                # Combina ambos genes en una sola figura
) & 
  theme_bw() &                 # Estética limpia
  theme(
    legend.position = "none",   # Quitamos la leyenda para ahorrar espacio
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "italic", size = 14) # Nombres de genes en itálica
  )

# 4. Guardar la figura
ggsave(filename = paste(path.guardar, "VlnPlot_Endothelial_Markers.pdf", sep="/"), plot = vln_plot, width = 8,height = 5)

# ----------------------------------------
# Violin Plot of EC Markers 
# ----------------------------------------

# 2. Definir los genes de interés (formato correcto: Primera mayúscula para ratón)
genes_endoteliales <- c("Pecam1","Pik3ca")

# 3. Generar el VlnPlot
vln_plot <- VlnPlot(
  data, 
  features = genes_endoteliales, 
  group.by = "Layer_1", 
  cols = colors_use,       # Aplicamos tu paleta homogenizada
  pt.size = 0.1,                # Tamaño de los puntos (células) sobre el violín
  combine = TRUE                # Combina ambos genes en una sola figura
) & 
  theme_bw() &                 # Estética limpia
  theme(
    legend.position = "none",   # Quitamos la leyenda para ahorrar espacio
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "italic", size = 14) # Nombres de genes en itálica
  )

# 4. Guardar la figura
ggsave(filename = paste(path.guardar, "VlnPlot_pecam_Pi3k_Markers.pdf", sep="/"), plot = vln_plot, width = 8,height = 5)


# ----------------------------------------
# DimPlot eGFP positive and negative manual 
# ----------------------------------------

metadata<-data@meta.data
umaps_dim<-data@reductions$umap@cell.embeddings
metadata<-cbind(metadata,umaps_dim)

# | Control 
ct<-TRUE
pheno<-"WT"
metadata_plot<-metadata[which(metadata$Phenotype == pheno),]
metadata_plot$CellOfInterest<-"0"
metadata_plot$CellOfInterest[which(metadata_plot$eGFP_pos==ct)] <- "1"

control <- ggplot(metadata_plot, aes(x = umap_1, y = umap_2,col=CellOfInterest,fill=CellOfInterest)) +
   # Primero las "0" (detrás, más pequeñas)
    geom_jitter(
        data = metadata_plot[metadata_plot$CellOfInterest == "0",],
        aes(x = umap_1, y = umap_2),
        shape = 21,
        size = 0.5,
        fill = "#EDEDED",
        color = "#EDEDED"
    ) +
    # Luego las "1" (delante, más grandes)
    geom_jitter(
        data = metadata_plot[metadata_plot$CellOfInterest == "1",],
        aes(x = umap_1, y = umap_2),
        shape = 21,
        size = 1,
        fill = "#73C098",
        color = "#059963"
    ) +
    theme_bw() + 
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none"
    ) +
    xlab("UMAP_1")+
    ylab("UMAP_2")+
    labs(title=pheno)+
  theme(plot.title = element_text(face = "bold",hjust = 0.5))

# Treated 
ct<-TRUE
pheno<-"MUT"
metadata_plot<-metadata[which(metadata$Phenotype == pheno),]
metadata_plot$CellOfInterest<-"0"
metadata_plot$CellOfInterest[which(metadata_plot$eGFP_pos==ct)] <- "1"

treated <- ggplot(metadata_plot, aes(x = umap_1, y = umap_2,col=CellOfInterest,fill=CellOfInterest)) +
    # Primero las "0" (detrás, más pequeñas)
    geom_jitter(
        data = metadata_plot[metadata_plot$CellOfInterest == "0",],
        aes(x = umap_1, y = umap_2),
        shape = 21,
        size = 0.5,
        fill = "#EDEDED",
        color = "#EDEDED"
    ) +
    # Luego las "1" (delante, más grandes)
    geom_jitter(
        data = metadata_plot[metadata_plot$CellOfInterest == "1",],
        aes(x = umap_1, y = umap_2),
        shape = 21,
        size = 1,
        fill = "#73C098",
        color = "#059963"
    ) +
    theme_bw() + 
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.position = "none"
    ) +
    xlab("UMAP_1")+
    ylab("UMAP_2")+
  labs(title=pheno)+
  theme(plot.title = element_text(face = "bold",hjust = 0.5))

control + treated
ggsave(paste(path.guardar,"DimPlot_eGFP_Values_Manual.pdf",sep="/"),width = 10,height = 5)

# ----------------------------------------
# Bar Plot of the eGFP proportions 
# ----------------------------------------

data<-SetIdent(data,value="eGFP_pos")
data_eGFP<-subset(data,ident="TRUE")

matchSCore2::summary_barplot(class.fac = data_eGFP$Layer_1,obs.fac =data_eGFP$Phenotype)+scale_fill_manual(values = colors_use)
ggsave(filename=paste(path.guardar,"Barplot_eGFP_WT_MUT.pdf",sep="/"),width=5,height=5)

tabla_frecuencias <- table(data_eGFP$Layer_1, data_eGFP$Phenotype)
# Si la quieres convertir en un data.frame estándar:
df_tabla <- as.data.frame(table(data_eGFP$Layer_1, data_eGFP$Phenotype))
openxlsx::write.xlsx(df_tabla,file.path(path.guardar,"DF_eGFP_Merge.xlsx"))


# ----------------------------------------
# Dot Plot Markers
# ----------------------------------------

markers<-readxl::read_excel(file.path(data_dir, "WT_MUT_markers.xlsx"),sheet="Hoja2")
#markers<-markers[1:32,1:2]
# Ordenar marcadores según cluster_order

markers_ordered <- markers %>%
  dplyr::mutate(
    excel_order = row_number(),
    CellType = factor(CellType, levels = cluster_order)
  ) %>%
  dplyr::arrange(CellType, excel_order)

genes_ordered <- unique(markers_ordered$Gene)


# --------------------------------------------------
# | DotPlot
# --------------------------------------------------

# Dotplot with RUIs version 

Bestholtz_palette <- c("#DEDAD6","#FEE392","#FEC44E","#FE9929","#ED6F11","#CC4C17","#993411","#65260C")

markers_ordered <- markers %>%
  dplyr::mutate(
    Gene = trimws(as.character(Gene)),
    CellType = trimws(as.character(CellType)),
    excel_order = dplyr::row_number(),
    CellType = factor(CellType, levels = cluster_order)
  ) %>%
  dplyr::arrange(CellType, excel_order)

gene_to_cluster <- markers_ordered %>%
  dplyr::select(
    gene = Gene,
    CellType,
    excel_order
  ) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(CellType, excel_order)

genes_ordered <- gene_to_cluster$gene

p <- DotPlot(
  data,
  group.by = "Layer_1",
  features = genes_ordered,
  col.min = 0,
  dot.scale = 5,
  scale = TRUE
)

p$data$features.plot <- factor(
  p$data$features.plot,
  levels = genes_ordered
)

p <- p +
  scale_x_discrete(
    limits = genes_ordered,
    drop = FALSE
  ) +
  theme(
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  scale_colour_gradientn(
    colours = Bestholtz_palette
  )

df <- gene_to_cluster %>%
  dplyr::transmute(
    x = factor(gene, levels = genes_ordered),
    CellType = as.character(CellType),
    y = 0.4,
    z = factor(CellType, levels = cluster_order)
  )

p_final <- p +
  geom_tile(
    data = df,
    aes(
      x = x,
      y = y,
      fill = z
    ),
    inherit.aes = FALSE,
    height = 0.3,
    width = 1,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = colors_use,
    breaks = cluster_order,
    drop = FALSE
  )

ggsave(plot=p_final,file=paste(path.guardar,"Rui_DotPlot_IdentityMarkers_Test_WT.pdf",sep="/"),width=20,height=5)
