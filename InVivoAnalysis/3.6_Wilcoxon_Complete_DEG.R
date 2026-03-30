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
library(clusterProfiler)
library(writexl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Wilcoxon/eGFP_WT_MUT_Total",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

# --- Carga y preparación ---
print("Load data")
# Path obj: 
path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")
data <- readRDS(paste(path.data,"Data_Subcluster.rds",sep="/"))
data$Phenotype<-factor(data$Phenotype,levels=c("WT","MUT"),labels=c("WT","MUT"))
DefaultAssay(data)<-"RNA"
cluster<-"0"
ident<-paste("SubCluster",cluster,"02",sep="_")
data<-SetIdent(data,value=ident)

# ----------------------------------------
# Annotation Cell Types 
# ----------------------------------------

data@meta.data$AnnotLayer<- revalue(data@meta.data[[ident]], c("0_0" = "C1",
                                                                      "0_1" = "C2",
                                                                      "1" = "C3",
                                                                      "2" = "C4",
                                                                      "3" = "C2",
                                                                      "4" = "C5"))
data<-SetIdent(data,value="AnnotLayer")
ident<-"AnnotLayer"
data$AnnotLayer<-factor(data$AnnotLayer,levels=c("C1","C2","C3","C4","C5")) 

data<-SetIdent(data,value="AnnotLayer")
universe<-rownames(data)
print("data loaded")

DefaultAssay(data) <- "RNA"

# Asegurar niveles
data$Phenotype <- factor(data$Phenotype, levels=c("WT","MUT"), labels=c("WT","MUT"))

# Split eGFP SOLO para MUT
data$eGFP_pos <- FetchData(data, vars="eGFP")$eGFP > 1.5
data$eGFP_pos <- factor(data$eGFP_pos, levels=c(FALSE, TRUE), labels=c("FALSE","TRUE"))
data$eGFP_pheno<-paste(data$Phenotype,data$eGFP_pos,sep="_")

ident<-"eGFP_pheno"
data<-SetIdent(data,value=ident)
data <- SetIdent(data, value = "eGFP_pheno")
markers <- FindMarkers(data,ident.1 = "MUT_TRUE",ident.2 = "WT_TRUE")
markers$gene <- rownames(markers)
markers$Classification <- ifelse(markers$avg_log2FC < 0, "Down", "Up")
writexl::write_xlsx(markers,paste(path.guardar,"Markers.xlsx",sep="/"))
saveRDS(markers,paste(path.guardar,"Markers.rds",sep="/"))

# | Estimate GO Terms 
markers_go<-markers[which(markers$Classification=="Up"),]
markers_go<-markers_go[which(markers_go$p_val_adj<0.05),]
UP <- enrichGO(gene=markers_go$gene,
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
writexl::write_xlsx(go_terms,paste(path.guardar,"GO_Terms.xlsx",sep="/"))
saveRDS(go_terms,paste(path.guardar,"GOTerms.rds",sep="/"))

# Enrichment Hallmark (mouse) -------------------------------------------------------------------------------------------------
hallmark_gene_sets_mm <- msigdbr(species = "Mus musculus", category = "H")

enrichment_results_hallmark <- enricher(
  gene = markers_go$gene,
  TERM2GENE = hallmark_gene_sets_mm[, c("gs_name", "gene_symbol")]
)

Hallmark <- enrichment_results_hallmark@result
writexl::write_xlsx(Hallmark, file.path(path.guardar, "Hallmark_mouse.xlsx"))
saveRDS(Hallmark, file.path(path.guardar, "Hallmark_mouse.rds"))


# Enrichment KEGG (mouse) -----------------------------------------------------------------------------------------------------
kegg_gene_sets_mm <- msigdbr(
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:KEGG"
)

enrichment_results_kegg <- enricher(
  gene = markers_go$gene,
  TERM2GENE = kegg_gene_sets_mm[, c("gs_name", "gene_symbol")]
)

KEGG <- enrichment_results_kegg@result

writexl::write_xlsx(KEGG, file.path(path.guardar, "KEGG_mouse.xlsx"))
saveRDS(KEGG, file.path(path.guardar, "KEGG_mouse.rds"))

# GSEA KEGG (mouse) -----------------------------------------------------------------------------------------------------

organism <-  "org.Mm.eg.db"

# # we want the log2 fold change 
original_gene_list <- markers$avg_log2FC
# # name the vector
names(original_gene_list) <- markers$gene

# ---- KEEG

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
 # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <- markers[markers$gene %in% dedup_ids$SYMBOL,]

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

openxlsx::write.xlsx(kk2@result,paste(path.guardar,"GSEA_KEEG.xlsx",sep="/"))
saveRDS(kk2@result, file.path(path.guardar, "GSEA_KEEG.rds"))

# ---- HALLMARK (MSigDB) GSEA (mouse)

# 1) Hallmark gene sets for mouse via msigdbr
#    (Requires: msigdbr, clusterProfiler, dplyr)
library(msigdbr)
library(clusterProfiler)
library(dplyr)

hallmark_mm <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

# TERM2GENE format required by clusterProfiler::GSEA
hallmark_term2gene <- hallmark_mm %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

# 2) Build ranked list (named numeric vector)
#    Assumes `markers` has columns: gene (SYMBOL) and avg_log2FC
hallmark_gene_list <- markers$avg_log2FC
names(hallmark_gene_list) <- markers$gene

# Remove NA + duplicated gene names (keep the one with largest absolute effect, or first)
hallmark_gene_list <- hallmark_gene_list[!is.na(hallmark_gene_list)]
hallmark_gene_list <- hallmark_gene_list[!duplicated(names(hallmark_gene_list))]

# Sort decreasing (required)
hallmark_gene_list <- sort(hallmark_gene_list, decreasing = TRUE)

# 3) Run GSEA against Hallmarks
hh <- GSEA(
  geneList       = hallmark_gene_list,
  TERM2GENE      = hallmark_term2gene,
  minGSSize      = 3,
  maxGSSize      = 800,
  pvalueCutoff   = 0.05,
  pAdjustMethod  = "none",   # keep as in your KEGG example; consider "BH" for FDR control
  verbose        = FALSE
)

# 4) Save results
openxlsx::write.xlsx(hh@result, file.path(path.guardar, "GSEA_HALLMARK.xlsx"))
saveRDS(hh@result, file.path(path.guardar, "GSEA_HALLMARK.rds"))