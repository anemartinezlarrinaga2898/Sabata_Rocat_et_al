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
library(readxl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.3_Dowstream_Analysis/Signatures",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

# Path obj: 
path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")
##################################################################################################################################

print("Load data")
data <- readRDS(paste(path.data,"Data_Anotado.rds",sep="/"))
data<-SetIdent(data,value="AnnotLayer")
print("data loaded")

# Load the excel with the signatures 
path.obj<-paste(path.guardar_original,"Figures/Summaries",sep="/")
signatures <- read_excel(file.path(path.obj,"Signatures_PI3K.xlsx"))
signatures$Gene<-stringr::str_to_title(signatures$Gene)
signature_list <- split(signatures,signatures$Signature)
gene_list <- lapply(signature_list, function(x) x$Gene)

data <- UCell::AddModuleScore_UCell(data, features = gene_list)
signature_names <- paste0(names(signature_list), "_UCell")

# | Total Cells 

VlnPlot(data,features = signature_names,sort = T,log=TRUE,raster=FALSE,group.by="AnnotLayer")
ggsave(filename = paste(path.guardar,"VLN_all_ct.png",sep="/"),width = 15,height = 15)

VlnPlot(data,features = signature_names,sort = T,log=TRUE,raster=FALSE,split.by="Phenotype",group.by="AnnotLayer")
ggsave(filename = paste(path.guardar,"VLN_all_ct_bycondition.png",sep="/"),width = 15,height = 15)

# Only eGFP+ cells 
data<-SetIdent(data,value="eGFP_Levels")
data_eGFP<-subset(data,ident="eGFP +")
data_eGFP<-SetIdent(data_eGFP,value="AnnotLayer")
VlnPlot(data_eGFP,features = signature_names,sort = T,log=TRUE,raster=FALSE,group.by="AnnotLayer")
ggsave(filename = paste(path.guardar,"VLN_eGFP.png",sep="/"),width = 15,height = 15)

VlnPlot(data_eGFP,features = signature_names,sort = T,log=TRUE,raster=FALSE,split.by="Phenotype",group.by="AnnotLayer")
ggsave(filename = paste(path.guardar,"VLN_eGFP_ByPheno.png",sep="/"),width = 15,height = 15)