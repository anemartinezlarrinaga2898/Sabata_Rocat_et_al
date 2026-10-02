################################################################################
# PUBLICATION CODE
# FILE: 5.0_Annotate_WT.R
# PURPOSE: Generate the final WT endothelial annotation and associated publication plots.
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
library(harmony)
library(Rcpp)
library(plyr)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- PATHS$results
path.guardar<-paste(path.guardar_original,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Final_Dist",sep="/")
dir.create(path.guardar,recursive=TRUE, showWarnings = FALSE)

# Load the script with the functions 
source_required("0.3_Util_SeuratPipeline.R")

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(PATHS$results,"0.4_WT_Analysis/EC_Subset/RemoveClusters/Res03_Clus67_Res05_239_Res05_Clus2/Test_PCA_Dimns/PCA_Dimns_24",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

# --------------------------------------------------------------------------------
# 1) Define eGFp Levels: 
# --------------------------------------------------------------------------------

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# --------------------------------------------------------------------------------
# 2) Annotate Clusters
# --------------------------------------------------------------------------------

ident<- "Harmony_Log_res.0.7"
data@meta.data$AnnotLayer<- revalue(data@meta.data[[ident]], c("0" = "Capillary Ecs",
                                                                      "1" = "Arterial-capillary Ecs",
                                                                      "2" = "Activated Angiogenic Ecs",
                                                                      "3" = "Angiogenic pre-arterial Ecs",
                                                                      "4" = "Proliferative 2",
                                                                      "5" = "Proliferative 1",
                                                                      "6" = "Arterial Ecs",
                                                                      "7" = "Capillary Ecs",
                                                                      "8" = "Venous Ecs",
                                                                      "9" = "Tip Ecs",
                                                                      "10" = "Remove",
                                                                      "11" = "Remove"))

data <- SetIdent(data,value="AnnotLayer")
data_final <-subset(data,ident="Remove",invert=TRUE)
saveRDS(data_final,file.path(path.guardar,"Data_WT_Anotado.rds"))

# ----------------------------------------
# UMAP All populations
cluster_order<-c("Tip Ecs","Angiogenic pre-arterial Ecs","Activated Angiogenic Ecs","Proliferative 1","Proliferative 2","Venous Ecs","Capillary Ecs","Arterial-capillary Ecs","Arterial Ecs")
colors_use<-c("#b37087","#F48FB1","#f4cedb","#b5dcfb","#90caf9","#c5e1a5","#D7B393","#FFE082","#FF8A65")
names(colors_use) <- cluster_order
# ----------------------------------------

data_final$AnnotLayer<-factor(data_final$AnnotLayer,levels=cluster_order)

# UMAP Total
DimPlot(data_final,reduction = "umap",group.by = "AnnotLayer",label = TRUE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_Label.pdf",sep="/"),width=7,height=7)
DimPlot(data_final,reduction = "umap",group.by = "AnnotLayer",label = FALSE, cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabel.pdf",sep="/"),width=7,height=7)
DimPlot(data_final,reduction = "umap",group.by = "AnnotLayer",label = FALSE, label.size = 5,cols =colors_use,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(paste(path.guardar,"UMAP_All_Cell_Types_NoLabelNoLegend.pdf",sep="/"),width=7,height=7)
# ----------------------------------------
# Violin Plot of EC Markers 
# ----------------------------------------

# 2. Definir los genes de interés (formato correcto: Primera mayúscula para ratón)
genes_endoteliales <- c("Pecam1", "Cdh5")

# 3. Generar el VlnPlot
vln_plot <- VlnPlot(
  data_final, 
  features = genes_endoteliales, 
  group.by = "AnnotLayer", 
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
# Expression
# ----------------------------------------

FeaturePlot(data_final, features = c("Pik3ca","Pik3r1","Pten"), min.cutoff = "q9", order = T,pt.size = 0.5,ncol=2)&NoAxes()
ggsave(filename = paste(path.guardar,"Fp_PI3K_Levels.pdf",sep="/"),width = 10,height = 10)

FeaturePlot(data_final, features = c("Cxcr4","Unc5b","Esm1","Nr2f2","Bmx"), min.cutoff = "q9", order = T,pt.size = 0.5,ncol=2)&NoAxes()
ggsave(filename = paste(path.guardar,"Fp_Markers_Levels.pdf",sep="/"),width = 10,height = 10)

VlnPlot(data_final,features = c("Pik3ca","Pik3r1","Pten"),sort = T,log=TRUE,raster=FALSE,group.by="AnnotLayer",ncol=3,col=colors_use)
ggsave(filename=paste(path.guardar,"Vln_PI3K_Levels.pdf",sep="/"),width=15,height = 7)

# ----------------------------------------
# eGFP Levels
# ----------------------------------------

# Density plot of the eGFP Levels ...............
egfp_expr <- FetchData(data_final, vars = "eGFP")
ggplot(egfp_expr, aes(x = eGFP)) +
  geom_density(fill = "#93c5fd", alpha = 0.6) +
  labs(
    title = "Density plot de la expresión de eGFP",
    x = "Expresión eGFP",
    y = "Densidad")
ggsave(paste(path.guardar,"Density_eGFP_General.pdf",sep="/"),width=5,height=5)

# # Density plot of the eGFP levels per cluster ......
df2 <- FetchData(data_final, vars = c("eGFP"))
df2$cluster <- Idents(data_final)

ggplot(df2, aes(x = eGFP, color = cluster)) +
  geom_density() +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Density plot de eGFP por cluster",
    x = "Expresión eGFP",
    y = "Densidad"
  )
ggsave(paste(path.guardar,"Density_eGFP_Cluster.pdf",sep="/"),width=5,height=5)

# ----------------------------------------
# Bar Plot of the eGFP proportions 
# ----------------------------------------

data_final<-SetIdent(data_final,value="eGFP_pos")
data_eGFP<-subset(data_final,ident="TRUE")

matchSCore2::summary_barplot(class.fac = data_eGFP$AnnotLayer,obs.fac =data_eGFP$Phenotype)+scale_fill_manual(values = colors_use)
ggsave(filename=paste(path.guardar,"Barplot_eGFP_WT_MUT.pdf",sep="/"),width=5,height=5)


tabla_frecuencias <- table(data_eGFP$AnnotLayer, data_eGFP$Phenotype)
# Si la quieres convertir en un data.frame estándar:
df_tabla <- as.data.frame(table(data_eGFP$AnnotLayer, data_eGFP$Phenotype))
openxlsx::write.xlsx(df_tabla,file.path(path.guardar,"DF_eGFP_WT.xlsx"))


DimPlot(data,reduction = "umap",group.by = "eGFP_Levels",label = FALSE, label.size = 5,cols =c("blue","GREY"),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = file.path(path.guardar,paste("DimPlot_Total_eGFP_Levels_",ident,".png",sep="_")),width = 6,height = 6)

data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# UMAP 
DimPlot <- DimPlot(data_final,group.by = "eGFP_Levels",pt.size = 0.5,raster=FALSE,cols=c("#73C098","#EDEDED")) &NoAxes()
ggsave(filename=paste(path.guardar,"Dimplot_eGFP_pos_Split.pdf",sep="/"),plot=DimPlot,width=5,height=5)

# --- Define it manually: 

metadata<-data_final@meta.data
umaps_dim<-data_final@reductions$umap@cell.embeddings
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

ggsave(filename=paste(path.guardar,"Dimplot_eGFP_pos_WT.pdf",sep="/"),plot=control,width=5,height=5)

# ------------------------------------------------------
# | Estimate TOP 50 genes per cluster 
# ------------------------------------------------------

# 1. Asegurar la identidad de los clusters
Idents(data) <- "AnnotLayer"

# 2. Encontrar todos los marcadores
all_markers <- FindAllMarkers(
  data, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# 3. Filtrar los Top 50 por cluster y mantenerlos en un solo dataframe
# Usamos el pipe para ordenar por cluster y luego por log2FC
top50_single_sheet <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  ungroup()

# # 4. Guardar en una única página de Excel
# # El archivo tendrá una sola pestaña llamada "Top50_Markers"
# write.xlsx(
#   top50_single_sheet, 
#   file = file.path(path.guardar, "Top50_Markers_AllClusters_SingleSheet.xlsx"),
#   sheetName = "Top50_Markers",
#   rowNames = FALSE
# )

# ----------------------------------------
# Dot Plot Markers
# ----------------------------------------

markers<-readxl::read_excel(PATHS$wt_final_annotation,sheet="Final")
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
  data_final,
  group.by = "AnnotLayer",
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

ggsave(plot=p_final,file=paste(path.guardar,"Rui_DotPlot_IdentityMarkers.pdf",sep="/"),width=20,height=5)
