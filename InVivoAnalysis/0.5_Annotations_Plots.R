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
library(readxl)
library(gridExtra)
getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.1_SeuratPipeline/Annotations",sep="/")
dir.create(path.guardar,recursive=TRUE)

# Load the script with the functions 
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

args = commandArgs(trailingOnly=TRUE)
################################################################################

path.obj<-paste(Lista_Paths_Main$Path_Guardar_Retinas,"0.1_SeuratPipeline/Harmony",sep="/")
data <- readRDS(paste(path.obj,"ID_LogNormalized.rds",sep="/"))
DefaultAssay(data)<-"RNA"
ident<-args[1]
data <- SetIdent(data,value =ident)

path.guardar<-paste(path.guardar_original,"0.1_SeuratPipeline/Annotations",ident,sep="/")
dir.create(path.guardar,recursive=TRUE)

# Annotate everything as signature .............................................

Signatures <- read_excel("/ijc/USERS/amartinezl/MGRAUPERA_18/Signatures.xlsx")
Signatures_List <- split(Signatures,Signatures$Signature)
names(Signatures_List) <- stringr::str_replace_all(names(Signatures_List)," ","_")
names(Signatures_List) <- stringr::str_to_title(names(Signatures_List))

Lista_Signature <- list()

for(i in seq_along(Signatures_List)){
  ct<-Signatures_List[[i]]
  m<-stringr::str_to_title(ct$Gene)
  Lista_Signature[[i]]<-m
  names(Lista_Signature)[i]<-names(Signatures_List)[i]
}

data <- UCell::AddModuleScore_UCell(data, features = Lista_Signature)
signature_names <- paste0(names(Lista_Signature), "_UCell")

Lista_FP<-list()
Lista_VL<-list()

for(i in seq_along(signature_names)){
  sig_names<-signature_names[i]
  print(sig_names)
  Fp<-FeaturePlot(data, features = sig_names, min.cutoff = "q9", order = T) &NoAxes()
  Vln<-VlnPlot(data, features =sig_names ,sort = T)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  Lista_FP[[i]]<-Fp
  Lista_VL[[i]]<-Vln
}

x<-do.call("grid.arrange",c(Lista_FP,ncol=5,nrow=2))
ggsave(plot=x,filename=paste(path.guardar,"FeaturePlot_Annotation.png",sep="/"),width=30,heigh=12)

x<-do.call("grid.arrange",c(Lista_VL,ncol=5,nrow=2))
ggsave(plot=x,filename=paste(path.guardar,"ViolinPlot_Annotation.png",sep="/"),width=30,heigh=12)

Lista_FP<-list()
Lista_VL<-list()

for(i in seq_along(signature_names)){
  sig_names<-signature_names[i]
  print(sig_names)
  file_name<-paste("Fp_",sig_names,".png",sep="")
  Fp<-FeaturePlot(data, features = sig_names, min.cutoff = "q9", order = T) &NoAxes()
  ggsave(plot=Fp,filename=paste(path.guardar,file_name,sep="/"),width=5,heigh=5)
  
  Vln<-VlnPlot(data, features =sig_names ,sort = T)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  file_name<-paste("Vln_",sig_names,".png",sep="")
  ggsave(plot=Vln,filename=paste(path.guardar,file_name,sep="/"),width=5,heigh=5)
}

# Division estimation 
s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <-  stringr::str_to_title(cc.genes$g2m.genes)

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,pt.size = 1,raster=FALSE,split.by = "Phase")&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype_SplitByPhase.png",sep="/"),width = 15,height = 7)

col <-  c("#95BDED","#F5A040","#C8C2F2")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype.png",sep="/"),width = 7,height = 7)

DimPlot(data,reduction = "umap",group.by = "Phase",label = TRUE, label.size = 5,pt.size = 0.5,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phase.png",sep="/"),width = 7,height = 7)

FeaturePlot(data, features = "S.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score.png",sep="/"),width = 7,height = 7)
FeaturePlot(data, features = "G2M.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score.png",sep="/"),width = 7,height = 7)

FeaturePlot(data, features = "S.Score", min.cutoff = "q9", order = T,split.by = "Phenotype") &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score_SplitByPhase.png",sep="/"),width = 10,height = 7)
FeaturePlot(data, features = "G2M.Score", min.cutoff = "q9", order = T,split.by = "Phenotype") &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score_SplitByPhase.png",sep="/"),width = 10,height = 7)

