################################################################################
# SCRIPT: Downsampling sensitivity analysis for WT versus MUT differential expression
# AUTHOR: Ane Martinez Larrinaga
#
# DESCRIPTION:
# Balances WT and MUT eGFP-positive cell numbers within a selected endothelial population and repeats Seurat differential-expression analysis across 1,000 random downsampling rounds.
#
# INPUT:
#   results/0.6_Join_MUT_WT/EC_Subset/Annotations/Data_EC_Annotated.rds
#
# OUTPUT:
#   Per-round marker results and combined downsampling results for the selected population
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
library(optparse)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.6_Join_MUT_WT/DownStream/DEG/Downsampling",sep="/")
dir.create(path.guardar, recursive = TRUE, showWarnings = FALSE)

# Load project utility functions 
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

# Path obj: 

# Script parameters
option_list <- list(make_option(c("-i", "--index"), type = "character", help = "ct indes"))
opt <- parse_args(OptionParser(option_list = option_list))
idx <- opt$i
if (is.null(idx) || !nzchar(idx) || is.na(suppressWarnings(as.integer(idx)))) {
  stop("Provide a valid 1-based cluster index with --index.")
}
levels_needed <- c("WT", "MUT") 
phen_col<- "Phenotype"
fc_col <- "avg_log2FC"     # en Seurat suele ser avg_log2FC (comprueba tu objeto)
padj_col <- "p_val_adj"    # en Seurat es p_val_adj
alpha <- 0.05
##################################################################################################################################

# Data loading and preparation
print("Load data")
# Path obj: 
path.guardar_original <- results_root
path.obj<-paste(path.guardar_original,"0.6_Join_MUT_WT/EC_Subset/Annotations",sep="/")
data <- readRDS(paste(path.obj,"Data_EC_Annotated.rds",sep="/"))
data<-SetIdent(data,value="Layer_1")
universe<-rownames(data)
print("data loaded")
DefaultAssay(data) <- "RNA"

# Define eGFP levels:
data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("TRUE","FALSE"),labels=c("TRUE","FALSE"))
data$eGFP_pos_num <- ifelse(data$eGFP_pos == "TRUE", 1, 0)
data@meta.data <- data@meta.data %>% mutate(eGFP_Levels = case_when(eGFP_pos == TRUE ~ "eGFP +",eGFP_pos == FALSE ~ "eGFP -"))
data$eGFP_Levels<-factor(data$eGFP_Levels , levels=c("eGFP +","eGFP -"),labels=c("eGFP +","eGFP -"))

# eGFP+ (filtrar por meta.data mejor que por ident)
data_positive <- subset(data, subset = eGFP_pos == TRUE)

# ahora sí: cluster basado en Layer_1
data_positive <- SetIdent(data_positive, value="Layer_1")
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

# Repeat the balanced comparison across 1,000 downsampling rounds
round<-1000
Lista_Results_FindMarkers<-list()
Round_Summary <- vector("list", length = round)

for(i in seq(1,round,1)){
    print(paste("Round",i,sep=" "))
    cells_wt_ds  <- sample(cells_wt,  n_star)
    cells_mut_ds <- sample(cells_mut, n_star)
    cells_keep <- c(cells_wt_ds, cells_mut_ds)

    # Subset balanced cells
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
