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
library(clusterProfiler)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03_Clus_7_8_Res03_Clus_3_5_Res_03_Clus6/Downstream/DEG/eGFP_Mut_WT",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

#args = commandArgs(trailingOnly=TRUE)
################################################################################

path.data<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"Figures",sep="/")
data <- readRDS(paste(path.data,"Data_Anotado.rds",sep="/"))
data$eGFP_pos <- FetchData(data, vars = "eGFP")$eGFP > 2
data$eGFP_pos <- factor(data$eGFP_pos , levels=c("FALSE","TRUE"),labels=c("FALSE","TRUE"))
data$eGFP_pheno<-paste(data$Phenotype,data$eGFP_pos,sep="_")
universe<-rownames(data)
DefaultAssay(data)<-"RNA"
data$Phenotype<-factor(data$Phenotype,levels=c("WT","MUT"),labels=c("WT","MUT"))

#ident<-args[1]
ident<-"AnnotLayer"
data<-SetIdent(data,value=ident)

# data<-SetIdent(data,value=ident)
# cluster<-"0"
# cluster_res_02_name<-paste("SubCluster",cluster,"02",sep="_")
# data<-FindSubCluster(data,cluster,"Harmony_Log",subcluster.name = cluster_res_02_name,resolution = 0.2,algorithm = 1)
# data<-SetIdent(data,value=cluster_res_02_name)

cluster<-unique(data@active.ident)

Lista_Markers<-list()

for(i in seq_along(cluster)){

    c <- cluster[i]

    data_c <- subset(data, ident = c)
    data_c <- SetIdent(data_c, value = "eGFP_pheno")

    # Contar células por grupo
    cell_counts <- table(Idents(data_c))

    # Comprobar que ambos grupos existen y tienen ≥ 3 células
    if(
        all(c("MUT_TRUE", "WT_TRUE") %in% names(cell_counts)) &&
        cell_counts["MUT_TRUE"] >= 3 &&
        cell_counts["WT_TRUE"] >= 3
    ){

        markers <- FindMarkers(
            data_c,
            ident.1 = "MUT_TRUE",
            ident.2 = "WT_TRUE"
        )

        markers$gene <- rownames(markers)
        markers$cluster <- c
        markers$Classification <- ifelse(markers$avg_log2FC < 0, "Down", "Up")

        Lista_Markers[[length(Lista_Markers) + 1]] <- markers
        names(Lista_Markers)[length(Lista_Markers)] <- paste("Cluster", c, sep = "_")

    } else {

        message(
            paste(
                "Skipping cluster", c,
                "- insufficient cells:",
                paste(names(cell_counts), cell_counts, collapse = ", ")
            )
        )
    }
}


file_name<-paste("DEG_Mut_WT_eGFP_Pos_",ident,".xlsx",sep="")
openxlsx::write.xlsx(Lista_Markers,paste(path.guardar,file_name,sep="/"))
