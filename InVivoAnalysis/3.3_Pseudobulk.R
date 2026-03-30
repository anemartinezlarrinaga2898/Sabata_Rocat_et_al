# SCRIPT: Harmony with LogWT_TRUEized data
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
library(clusterProfiler)
library(optparse)
library(DESeq2)
library(writexl)
library(openxlsx)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

#========================
# 0) Parámetros
#========================
assay_use <- "RNA"
min_cells <- 30
min_pct   <- 0.01  # 1%

option_list <- list(make_option(c("-c", "--celltype"), type = "character", help = "Number of the tissue to process"))
opt <- parse_args(OptionParser(option_list = option_list))

# Nombres de columnas en meta.data (ajusta si hace falta)
col_patient <- "ID"
col_condition <- "Phenotype"        # WT_TRUE / MUT_TRUE
col_doublet <- "DoubletClassification_DoubletFinder"     # Singlet / Doublet
col_layer   <- "AnnotLayer"
cell_type <- opt$c

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/DEG/Pseudobulk/eGFP_WT_MUT",cell_type,sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))
################################################################################

# --- Carga y preparación ---
print("Load data")
# Path obj: 
path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")
data <- readRDS(paste(path.data,"Data_Anotado.rds",sep="/"))
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

ident<-"AnnotLayer"
data<-SetIdent(data,value=ident)

#========================
# 1) Filtrar singlets
#========================
DefaultAssay(data) <- assay_use

data <- SetIdent(data, value = col_doublet)
data<-subset(data,idents=c("Singlet"))

# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

data$samples <- paste(data$eGFP_pheno, data$ID,sep="_")
data$cell<-data[[col_layer]]

cts <- AggregateExpression(data, 
                    group.by = c("cell", "samples"),
                    assays = 'RNA',
                    slot = "counts",
                    return.seurat = FALSE)
# transpose
cts<-cts$RNA
cts.t <- t(as.matrix(cts))
# convert to data.frame
cts.t <- as.data.frame(cts.t)
# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))
# split data.frame
cts.split <- split.data.frame(cts.t,
                 f = factor(splitRows))
# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

# Let's run DE analysis 
# 1. Get counts matrix

idx<-which(names(cts.split.modified)==cell_type)

counts_ct <- cts.split.modified[[idx]]
storage.mode(counts_ct) <- "integer"

# 2. generate sample level metadata
colData <- data.frame(samples = colnames(counts_ct))

colData <- colData %>%
  filter(grepl("TRUE", samples)) %>%   # solo TRUE
  mutate(condition = case_when(
    grepl("^MUT-TRUE", samples) ~ "MUT",
    grepl("^WT-TRUE",  samples) ~ "WT")) %>%
  column_to_rownames(var = "samples")

colData$condition <- relevel(factor(colData$condition), ref = "WT")
colData$samples <- rownames(colData)

samples_keep <- colData$samples   # antes de hacer column_to_rownames()
count_matrix <- counts_ct[, samples_keep]

# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                       colData = colData,
                       design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)

# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res <- results(dds, name = "condition_MUT_vs_WT")
res$gene<-rownames(res)
saveRDS(res,paste(path.guardar,"Res_Repro.rds",sep="/"))
write.xlsx(res,file = paste(path.guardar,"Res_Repro.xlsx",sep="/"),rowNames = FALSE)

#========================
# 6) Enrichment Analysis
#========================

library(clusterProfiler)
library(AnnotationDbi)
library(msigdbr)

markers_filtered  <- res[res$log2FoldChange > 0,]
markers_filtered  <- markers_filtered[which(markers_filtered$padj<0.05),]

go_terms <- enrichGO(gene=markers_filtered$gene,
                    universe = universe,
                    OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE,
                    pAdjustMethod = "BH",
                    ont="BP",
                    minGSSize=5)
go_terms<-go_terms@result
writexl::write_xlsx(go_terms,paste(path.guardar,"GoTerms_Deseq2_Repro_SeuratFunction.xlsx",sep="/"))

# Perfrom also GSEA analysis 

# | Perfrom a GSEA anlaysis as we dont find any sig genes 

geneList <- res$stat
names(geneList) <- rownames(res)
geneList <- sort(geneList, decreasing = TRUE)

ego_gsea <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Mm.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 30,
  maxGSSize    = 200,
  pvalueCutoff = 0.1
)
saveRDS(ego_gsea,paste(path.guardar,"ego_gsea.rds",sep="/"))
write.xlsx(ego_gsea@result,file = paste(path.guardar,"GSEA_GO_results.xlsx",sep="/"),rowNames = FALSE)

res_sel <- ego_gsea@result |>
  dplyr::filter(p.adjust < 0.1) |>
  dplyr::filter(abs(NES) >= 1.3) |>
  dplyr::filter(setSize >= 20, setSize <= 300) |>
  dplyr::arrange(dplyr::desc(abs(NES)))
top15 <- head(res_sel, 15)
openxlsx::write.xlsx(top15, "GSEA_selected_top15.xlsx", rowNames = FALSE)