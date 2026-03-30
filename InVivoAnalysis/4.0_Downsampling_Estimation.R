# SCRIPT: Downsampling analyis
# AUTOR: ANE 
# Date: 16-01-2026

# Queremos ver si realmente el numero de DEG de las arterias se ve afectado por el numero de celulas, es decir si todos lo  hacemos downgrade y lo igualamos
# obtenemos el mismo patron? 

# Luego con el pseudobulk confirmaremos tambien que realmente que no pasa nada por comparar 80 vs 1000

##################################################################################################################################
set.seed(123)
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

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Downsampling",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

# Path obj: 
path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")

# | Parámetros del script
option_list <- list(make_option(c("-i", "--index"), type = "character", help = "ct indes"))
opt <- parse_args(OptionParser(option_list = option_list))
idx <- opt$i
levels_needed <- c("WT", "MUT") 
phen_col<- "Phenotype"
fc_col <- "avg_log2FC"     # en Seurat suele ser avg_log2FC (comprueba tu objeto)
padj_col <- "p_val_adj"    # en Seurat es p_val_adj
alpha <- 0.05
##################################################################################################################################
print("Load data")
data <- readRDS(paste(path.data,"Data_Anotado.rds",sep="/"))
data<-SetIdent(data,value="AnnotLayer")
print("data loaded")

# Define eGFP levels:
data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# eGFP+ (filtrar por meta.data mejor que por ident)
data_positive <- subset(data, subset = eGFP_pos == TRUE)

# ahora sí: cluster basado en AnnotLayer
data_positive <- SetIdent(data_positive, value="AnnotLayer")
cluster <- levels(Idents(data_positive))  # o sort(unique(Idents(data_positive)))
cluster_to_study <- as.character(cluster[as.integer(idx)])
print(paste("Markers from:",cluster_to_study,sep=" "))
data_c <- subset(data_positive, ident = cluster_to_study)

md <- data_c@meta.data
cells_wt  <- rownames(md)[as.character(md[[phen_col]]) == levels_needed[1]]
cells_mut <- rownames(md)[as.character(md[[phen_col]]) == levels_needed[2]]

n_wt  <- length(cells_wt)
n_mut <- length(cells_mut)

n_star <- min(n_wt, n_mut)   # tamaño balanceado
if (n_star < 10) stop("Muy pocas células para downsampling robusto en este cluster")

# Foor loop to perfom 1000 rounds
round<-1000
Lista_Results_FindMarkers<-list()
Round_Summary <- vector("list", length = round)

for(i in seq(1,round,1)){
    print(paste("Round",i,sep=" "))
    cells_wt_ds  <- sample(cells_wt,  n_star)
    cells_mut_ds <- sample(cells_mut, n_star)
    cells_keep <- c(cells_wt_ds, cells_mut_ds)

    # Subse4t the cells
    data_c_pos_balanced <- subset(data_c, cells = cells_keep)
    data_c_pos_balanced<-SetIdent(data_c_pos_balanced,value="Phenotype")
    markers_true <- FindMarkers(data_c_pos_balanced,ident.1 = "MUT",ident.2 = "WT",logfc.threshold = 0,min.pct = 0)
    markers_true$Round<-paste("Round_",i,sep="")
    Lista_Results_FindMarkers[[i]]<-markers_true

    # ---- resumen de conteos por ronda ----
  df <- markers_true %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(
      significant = .data[[padj_col]] < alpha,
      direction = case_when(
        .data[[fc_col]] > 0 ~ "up",
        .data[[fc_col]] < 0 ~ "down",
        TRUE ~ "zero"
      )
    )

  Round_Summary[[i]] <- tibble::tibble(
    Round = paste0("Round_", i),
    n_total = nrow(df),
    n_significant = sum(df$significant, na.rm = TRUE),
    n_not_significant = sum(!df$significant, na.rm = TRUE),
    n_up = sum(df$direction == "up", na.rm = TRUE),
    n_down = sum(df$direction == "down", na.rm = TRUE),
    n_up_significant = sum(df$direction == "up" & df$significant, na.rm = TRUE),
    n_down_significant = sum(df$direction == "down" & df$significant, na.rm = TRUE),
    n_up_not_significant = sum(df$direction == "up" & !df$significant, na.rm = TRUE),
    n_down_not_significant = sum(df$direction == "down" & !df$significant, na.rm = TRUE)
  )
}

# Save Results
file_name<-paste("Lista_Markers_Repro_",cluster_to_study,".rds",sep="")
saveRDS(Lista_Results_FindMarkers,file.path(path.guardar,file_name))

# Merge results in a single df
DF_complete<-do.call(rbind,Lista_Results_FindMarkers)
DF_complete$Cluster<-cluster_to_study
file_name<-paste("DF_Complete_",cluster_to_study,".rds",sep="")
saveRDS(DF_complete,file.path(path.guardar,file_name))
