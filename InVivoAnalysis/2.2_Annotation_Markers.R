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

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.2_EC_Subset/Subset_Clus_Res_03/GeneralAnnot_Plots",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.2_EC_Subset/Subset_Clus_Res_03/Harmony",sep="/")
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
DefaultAssay(data)<-"RNA"

# | Markers Definition: 

endothelial_markers <- list(
  vein_general = c("Aqp1","Adgrg6","Mafb","Itih5","Slco1a4","Olfml2a","Col23a1","Nox4"),
  cap_general = c("Mfsd2a","Cxcl12","Gja4","Glul","Slc16a1","Pdgfd","Slco1c1","Apcdd1","Slc38a5"),
  proliferating_general = c("Top2a","Mki67","Ccnb1","Cdk1","Cenpf","Prc1"),
  tip_general = c("Apln","Ednrb","Angpt2","Tnfrsf9","Igfbp3","Lgals1"),
  artery_general = c("Efnb2","Gja5","Unc5b","Hey1","Sox17","Fbln5","Tmem100"),
  vein_lymphatic=c("Nos3","Nr2f2","Gpr182","Lama1","Npnt","Flvcr2","Angptl4","Ramp2","Pltp")
)

# | Feature Plot of the genes 

for(i in seq_along(endothelial_markers)){
    features<-endothelial_markers[[i]]
    names_features<-names(endothelial_markers)[i]

    FeaturePlot(data, features = features, min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
    file_name<-paste("Fp_",names_features,".png",sep="")
    ggsave(filename = paste(path.guardar,file_name,sep="/"),width = 10,height = 10)
}

# | Dot Plot of the genes 

# Res 0.5
data <- readRDS(paste(path.obj,"Harmony.rds",sep="/"))
data<-SetIdent(data,value="Harmony_Log_res.0.5")

for(i in seq_along(endothelial_markers)){

    features<-endothelial_markers[[i]]
    names_features<-names(endothelial_markers)[i]

    DotPlot(data, features = endothelial_markers)&
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))
    file_name<-paste("DotPlot_",names_features,".png",sep="")
    ggsave(filename = paste(path.guardar,file_name,sep="/"),width = 10,height = 10)
}

# DotPlot all together 

DotPlot(data, features = endothelial_markers)&
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=rel(1.1), face = "plain"),
        axis.text.y = element_text(size = rel(1.05), face = "plain"))
ggsave(filename = paste(path.guardar,"DotPlot_Total.png",sep="/"),width = 15,height = 7)
