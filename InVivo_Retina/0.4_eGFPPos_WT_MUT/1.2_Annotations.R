################################################################################
# SCRIPT: Annotation of eGFP-positive endothelial cell populations
# AUTHOR: Ane Martinez Larrinaga
# DESCRIPTION:
#   Assigns endothelial population labels from Harmony clusters and generates
#   annotated UMAP and marker DotPlot visualizations.
# INPUT:
#   0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony/Harmony.rds
# OUTPUT:
#   Annotation figures under 0.6_Join_MUT_WT/EC_Subset/eGFP/FinalAnnotations/.
# NOTE:
#   The original script did not save the annotated Seurat object because the
#   saveRDS() line was commented out. This behavior is preserved.
################################################################################

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), mustWork = FALSE)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(plyr)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- project_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/eGFP/FinalAnnotations",sep="/")
dir.create(path.guardar,recursive = TRUE, showWarnings = FALSE)

# Load the script with the functions 
source(file.path(project_root, "utils", "0.3_Util_SeuratPipeline.R"))

################################################################################

path.obj<-paste(project_root,"0.6_Join_MUT_WT/EC_Subset/eGFP/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

ident<-"Harmony_Log_res.0.7"
data@meta.data$Layer_1<- revalue(data@meta.data[[ident]], c("0" = "Proliferative 1",
                                                               "1" = "VenousCapillary EC",
                                                               "2" = "VenousCapillary EC",
                                                               "3" = "Proliferative 2",
                                                               "4" = "Proliferative 1",
                                                               "5" = "ArteryCapillary EC",
                                                               "6" = "Angiogenic EC",
                                                               "7" = "Capillary",
                                                               "8" = "VenousCapillary EC",
                                                               "9" = "Arterial EC"))

# Define colors: 
# ----------------------------------------
# UMAP All populations
cluster_order<-c("Angiogenic EC","Proliferative 1","Proliferative 2","VenousCapillary EC","Capillary","ArteryCapillary EC","Arterial EC")
colors_use<-c("#f4cedb","#b5dcfb","#90caf9","#C5E1A5","#e1be8f","#FFE082","#FF8A65")
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

df <- data.frame(
  Gene = c(
    "Esm1", "Angpt2", "Kcne3", "Odc1", "Apln", 
    "Ccne2", "Top2a",
    "Cenpa", "Ccnb2", 
    "Adgrg6", "Pcdh7", "Nr2f2", 
    "Hmcn1", "Dach1", "Mfsd2a", 
    "Apod","Gja4",  "Unc5b", "Cxcr4", "Hey1", 
    "Mgp", "Bmx", "Gja5"
  ),
  CellType = c(
    "Angiogenic EC", "Angiogenic EC", "Angiogenic EC", "Angiogenic EC", "Angiogenic EC", 
    "Proliferative 1", "Proliferative 1", 
    "Proliferative 2", "Proliferative 2", 
    "VenousCapillary EC", "VenousCapillary EC", "VenousCapillary EC", 
    "Capillary", "Capillary", "Capillary", 
    "ArteryCapillary EC", "ArteryCapillary EC", "ArteryCapillary EC", "ArteryCapillary EC", "ArteryCapillary EC", 
    "Arterial EC", "Arterial EC", "Arterial EC"
  ),
  stringsAsFactors = FALSE
)

# Dotplot with RUIs version 

Bestholtz_palette <- c("#DEDAD6","#FEE392","#FEC44E","#FE9929","#ED6F11","#CC4C17","#993411","#65260C")

markers_ordered <- df %>%
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

ggsave(plot=p_final,file=paste(path.guardar,"Rui_DotPlot_IdentityMarkers.pdf",sep="/"),width=15,height=5)
#saveRDS(data,file.path(path.guardar,"Data_Anotado.rds"))