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
library(plyr)
library(clusterProfiler)
library(msigdbr)
library(writexl)
library(optparse)
library(AnnotationDbi)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Wilcoxon/GoTerms_PerCluster",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

option_list <- list(make_option(c("-i", "--index"), type = "integer", help = "Number of the tissue to process"))
opt <- parse_args(OptionParser(option_list = option_list))
idx<-opt$i
################################################################################

# --- Carga y preparación ---
print("Load data")
# Path obj: 
path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")
data <- readRDS(paste(path.data,"Data_Subcluster.rds",sep="/"))
universe<-rownames(data)

path_file<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Wilcoxon/eGFP_WT_MUT",sep="/")
markers<-do.call(rbind,readRDS(file.path(path_file,"DEG_Mut_WT_eGFP_Pos_AnnotLayer.rds")))
cluster<-unique(markers$cluster)

c<-cluster[idx]
marker_c<-markers[which(markers$cluster==c),]
marker_c<-marker_c[which(marker_c$p_val_adj < 0.05),]
marker_c<-marker_c[which(marker_c$Classification == "Up"),]

# ----------------------------------------
# | Estimate GO Terms 
# ----------------------------------------

UP <- enrichGO(gene=marker_c$gene,
                    universe = universe,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE,
                    pAdjustMethod = "BH",
                    ont="BP",
                    minGSSize=10)

go_terms<-UP@result
file_name<-paste(c,"_GO_Terms.xlsx",sep="")
writexl::write_xlsx(go_terms,paste(path.guardar,file_name,sep="/"))
file_name<-paste(c,"_GO_Terms.rds",sep="")
saveRDS(go_terms,paste(path.guardar,file_name,sep="/"))

# ----------------------------------------
# | Estimate Hallmarks
# ----------------------------------------

hallmark_gene_sets_mm <- msigdbr(species = "Mus musculus", category = "H")

enrichment_results_hallmark <- enricher(
  gene = marker_c$gene,
  TERM2GENE = hallmark_gene_sets_mm[, c("gs_name", "gene_symbol")]
)

Hallmark <- enrichment_results_hallmark@result
file_name<-paste(c,"_Hallmark_mouse.xlsx",sep="")
writexl::write_xlsx(Hallmark, file.path(path.guardar, file_name))
file_name<-paste(c,"_Hallmark_mouse.rds",sep="")
saveRDS(Hallmark, file.path(path.guardar, file_name))

# ----------------------------------------
# | Estimate KEEG
# ----------------------------------------

kegg_gene_sets_mm <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:KEGG")

enrichment_results_kegg <- enricher(
  gene = marker_c$gene,
  TERM2GENE = kegg_gene_sets_mm[, c("gs_name", "gene_symbol")])

KEGG <- enrichment_results_kegg@result

file_name<-paste(c,"_KEGG_mouse.xlsx",sep="")
writexl::write_xlsx(KEGG, file.path(path.guardar, file_name))
file_name<-paste(c,"_KEGG_mouse.rds",sep="")
saveRDS(KEGG, file.path(path.guardar, file_name))

# ----------------------------------------
# | GSEA
# ----------------------------------------
organism <-  "org.Mm.eg.db"

# --- GO Terms
c<-cluster[idx]
marker_c<-markers[which(markers$cluster==c),]

# # we want the log2 fold change 
original_gene_list <- marker_c$avg_log2FC
# # name the vector
names(original_gene_list) <- marker_c$gene

gene_list<-na.omit(original_gene_list)
# # sort the list in decreasing order (required for clusterProfiler)
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

file_name<-paste(c,"_GSEA_GOTerms.xlsx",sep="")
openxlsx::write.xlsx(gse@result,paste(path.guardar,file_name,sep="/"))
file_name<-paste(c,"_GSEA_GOTerms.rds",sep="")
saveRDS(gse@result, file.path(path.guardar, file_name))

# ---- KEEG

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
 # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <- marker_c[marker_c$gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y <- dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

kegg_organism <- "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

file_name<-paste(c,"_GSEA_KEEG.xlsx",sep="")
openxlsx::write.xlsx(kk2@result,paste(path.guardar,file_name,sep="/"))
file_name<-paste(c,"_GSEA_KEEG.rds",sep="")
saveRDS(kk2@result, file.path(path.guardar, file_name))