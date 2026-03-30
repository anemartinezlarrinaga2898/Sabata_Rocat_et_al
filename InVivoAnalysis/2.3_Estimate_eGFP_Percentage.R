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
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony/eGFP_Percentage",sep="/")
#path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Harmony/eGFP_Percentage",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony",sep="/")
#path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Harmony",sep="/")
data <- readRDS(paste(path.data,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"
data$Phenotype<-factor(data$Phenotype,levels=c("WT","MUT"),labels=c("WT","MUT"))

ident<-args[1]
data<-SetIdent(data,value=ident)
cluster<-"0"
cluster_res_02_name<-paste("SubCluster",cluster,"02",sep="_")
data<-FindSubCluster(data,cluster,"Harmony_Log",subcluster.name = cluster_res_02_name,resolution = 0.2,algorithm = 1)
data<-SetIdent(data,value=cluster_res_02_name)
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Harmony/eGFP_Percentage",cluster_res_02_name,sep="/")
dir.create(path.guardar,recursive=TRUE)

# PI3K Levels

FeaturePlot(data, features = c("Pik3ca","Pik3r1"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(filename = paste(path.guardar,"Fp_PI3K_Levels.png",sep="/"),width = 10,height = 10)

VlnPlot(data,features = c("Pik3ca","Pik3r1"),sort = T,log=TRUE,raster=FALSE)
ggsave(filename=paste(path.guardar,"Vln_PI3K_Levels.png",sep="/"),width=15,heigh=7)

VlnPlot(data,features = c("Pik3ca","Pik3r1"),sort = T,log=TRUE,raster=FALSE,split.by="Phenotype")+theme(legend.position="bottom")
ggsave(filename=paste(path.guardar,"Vln_PI3K_Levels_Pheno.png",sep="/"),width=15,heigh=7)

# Dimplot split by phenotype 

DimPlot <- DimPlot(data,group.by = ident,pt.size = 1,raster=FALSE,split.by="Phenotype") &NoAxes()
ggsave(filename=paste(path.guardar,"Dimplot_ByPheno_Split.png",sep="/"),plot=DimPlot,width=10,height=5)

# | DensityPlot of eGFP Expression 

egfp_expr <- FetchData(data, vars = "eGFP")
ggplot(egfp_expr, aes(x = eGFP)) +
  geom_density(fill = "#93c5fd", alpha = 0.6) +
  labs(
    title = "Density plot de la expresión de eGFP",
    x = "Expresión eGFP",
    y = "Densidad")
ggsave(paste(path.guardar,"Density_eGFP_General.png",sep="/"),width=5,height=5)

# | DensityPlot por cluster 

df2 <- FetchData(data, vars = c("eGFP", "Phenotype"))
df2$cluster <- Idents(data)

ggplot(df2, aes(x = eGFP, color = cluster)) +
  geom_density() +
  facet_wrap(~ cluster, scales = "free_y") +
  labs(
    title = "Density plot de eGFP por cluster",
    x = "Expresión eGFP",
    y = "Densidad"
  )
ggsave(paste(path.guardar,"Density_eGFP_Cluster.png",sep="/"),width=5,height=5)

# Density Plot por phenotype 

df <- FetchData(data, vars = c("eGFP", "Phenotype"))

ggplot(df, aes(x = eGFP, fill = Phenotype, color = Phenotype)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Density plot de eGFP por condición",
    x = "Expresión eGFP",
    y = "Densidad"
  ) +
  scale_fill_manual(values = c("WT" = "#60a5fa", "MUT" = "#f87171")) +
  scale_color_manual(values = c("WT" = "#3b82f6", "MUT" = "#ef4444")) 
ggsave(paste(path.guardar,"Density_eGFP_Phenotype.png",sep="/"),width=5,height=5)


# Define the total eGFP cells 

data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 2

data$eGFP_pos <- factor(data$eGFP_pos , levels=c("FALSE","TRUE"),labels=c("FALSE","TRUE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)

DimPlot <- DimPlot(data,group.by = "eGFP_pos",pt.size = 0.5,raster=FALSE,split.by="Phenotype",cols=c("grey","blue")) &NoAxes()
ggsave(filename=paste(path.guardar,"Dimplot_eGFP_pos_Split.png",sep="/"),plot=DimPlot,width=10,height=5)

# Estimate the percentaje of eGFP cells per cluster 

df_clusters <- data@meta.data %>%
  dplyr::mutate(cluster = Idents(data)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    total = n(),
    eGFP_pos = sum(eGFP_pos_num),
    pct_eGFP = 100 * eGFP_pos / total
  )

# Estimate the percentaje of eGFP per phenotype 

df_conditions <- data@meta.data %>%
  group_by(Phenotype) %>%
  summarise(
    total = n(),
    eGFP_pos = sum(eGFP_pos_num),
    pct_eGFP = 100 * eGFP_pos_num / total
  )

# Estimate the percentaje of eGFP per cluster + phenotype  

df_cluster_cond <- data@meta.data %>%
  mutate(cluster = Idents(data)) %>%
  group_by(cluster, Phenotype) %>%
  summarise(
    total = n(),
    eGFP_pos = sum(eGFP_pos_num),
    pct_eGFP = 100 * eGFP_pos_num / total
  )

ggplot(df_cluster_cond, aes(x = Phenotype, y = cluster, fill = pct_eGFP)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#9BCD9B") +
  geom_text(aes(label = sprintf("%.1f%%", pct_eGFP)), color = "black", size = 3) +
  labs(
    title = "",
    x = "Condition",
    y = "Cluster",
    fill = "%"
  ) +
  theme(
    panel.background = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    legend.background = element_blank())
ggsave(paste(path.guardar,"Heatmap_Percentage_eGFP_Cluster_Phenotype.png",sep="/"),width=5,height=5)


# Extra plot: 

# ---------------------------------------------
# STACKED BARPLOT: eGFP+ vs eGFP- por cluster
# ---------------------------------------------

df_stack <- data@meta.data %>%
  mutate(cluster = Idents(data)) %>%
  group_by(cluster, eGFP_pos) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pct = 100 * n / sum(n))

ggplot(df_stack, aes(x = cluster, y = pct, fill = eGFP_pos)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "steelblue")) +
  labs(
    title = "Porcentaje de eGFP+ / eGFP− por cluster",
    x = "Cluster",
    y = "Porcentaje",
    fill = "eGFP+"
  ) 

ggsave(paste(path.guardar, "Stacked_eGFP_byCluster.png", sep="/"), width=7, height=5)

# ---------------------------------------------
# BARPLOT SIMPLE: % eGFP+ por cluster
# ---------------------------------------------

ggplot(df_clusters, aes(x = cluster, y = pct_eGFP)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", pct_eGFP)), vjust = -0.3) +
  labs(
    title = "% de células eGFP+ por cluster",
    x = "Cluster",
    y = "Porcentaje eGFP+"
  ) 

ggsave(paste(path.guardar, "Barplot_pct_eGFP_byCluster.png", sep="/"), width=7, height=5)

# -------------------------------------------------------------------------
# STACKED BARPLOT: eGFP+ vs eGFP– por cluster y por phenotype
# -------------------------------------------------------------------------

df_cluster_pheno <- data@meta.data %>%
  mutate(
    cluster = Idents(data),
    group = paste(Phenotype, eGFP_pos, sep = "_")
  ) %>%
  group_by(cluster, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(pct = 100 * n / sum(n))

ggplot(df_cluster_pheno, aes(x = cluster, y = pct, fill = group)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Composición por cluster combinando Phenotype + eGFP",
    x = "Cluster",
    y = "Porcentaje",
    fill = "Grupo"
  ) +
  scale_fill_manual(values = c(
    "WT_FALSE"  = "#a7c7e7",
    "WT_TRUE"   = "#3b82f6",
    "MUT_FALSE" = "#fca5a5",
    "MUT_TRUE"  = "#ef4444"
  )) 

ggsave(paste(path.guardar, "Stacked_Cluster_Phenotype_eGFP.png", sep="/"), width=8, height=5)

# | Save the tables 

wb <- createWorkbook()

addWorksheet(wb, "por_cluster")
writeData(wb, "por_cluster", df_clusters)

addWorksheet(wb, "por_condicion")
writeData(wb, "por_condicion", df_conditions)

addWorksheet(wb, "cluster_x_condicion")
writeData(wb, "cluster_x_condicion", df_cluster_cond)

saveWorkbook(wb, paste(path.guardar,"eGFP_results.xlsx",sep="/"), overwrite = TRUE)

# | Number of cells 

n_total <- nrow(data@meta.data)
df_n_total <- data.frame(total_cells = n_total)

df_n_condition <- data@meta.data %>%
  group_by(Phenotype) %>%
  summarise(total_cells = n())

df_n_cluster <- data@meta.data %>%
  mutate(cluster = Idents(data)) %>%
  group_by(cluster) %>%
  summarise(total_cells = n())

df_cells_cluster_condition <- data@meta.data %>%
  mutate(cluster = Idents(data)) %>%
  group_by(cluster, Phenotype) %>%
  summarise(total_cells = n())

wb2 <- createWorkbook()

addWorksheet(wb2, "total_global")
writeData(wb2, "total_global", df_n_total)

addWorksheet(wb2, "por_condicion")
writeData(wb2, "por_condicion", df_n_condition)

addWorksheet(wb2, "por_cluster")
writeData(wb2, "por_cluster", df_n_cluster)

addWorksheet(wb2, "cluster_x_condicion")
writeData(wb2, "cluster_x_condicion", df_cells_cluster_condition)

saveWorkbook(wb2, paste(path.guardar,"cell_counts.xlsx",sep="/"), overwrite = TRUE)