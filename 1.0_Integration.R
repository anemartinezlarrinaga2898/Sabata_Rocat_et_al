###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library(Seurat)
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)
library(gridExtra)
library(patchwork)
library(harmony)
library(Rcpp)

getPalette <-  colorRampPalette(brewer.pal(8, "Set1"))
col <-  getPalette(20)
path.guardar <- "0.2_SeuratPipeline/Integration"
dir.create(path.guardar)

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <-  stringr::str_to_title(cc.genes$g2m.genes)

Annotation_Function <- function(data,Lista_Markers,path.guardar){
  
  features_df<-do.call(rbind,Lista_Markers)
  cell_type<-names(Lista_Markers)
  
  Lista_Signature <- list()
  
  for(i in seq_along(Lista_Markers)){
    ct<-cell_type[i]
    print(ct)
    m<-features_df$Features[which(features_df$CellType==ct)]
    Lista_Signature[[i]]<-m
    names(Lista_Signature)[i]<-ct
  }
  
  data <- UCell::AddModuleScore_UCell(data, features = Lista_Signature)
  signature_names <- paste0(names(Lista_Signature), "_UCell")
  
  Lista_FP<-list()
  Lista_VL<-list()
  
  for(i in seq_along(signature_names)){
    sig_names<-signature_names[i]
    print(sig_names)
    Fp<-FeaturePlot(data, features = sig_names, min.cutoff = "q9", order = T) &NoAxes()
    Vln<-VlnPlot(data, features =sig_names ,sort = T,cols=col)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    Lista_FP[[i]]<-Fp
    Lista_VL[[i]]<-Vln
  }
  
  x<-do.call("grid.arrange",c(Lista_FP,ncol=5,nrow=2))
  ggsave(plot=x,filename=paste(path.guardar,"FeaturePlot_Annotation.png",sep="/"),width=30,heigh=12)
  
  x<-do.call("grid.arrange",c(Lista_VL,ncol=5,nrow=2))
  ggsave(plot=x,filename=paste(path.guardar,"ViolinPlot_Annotation.png",sep="/"),width=30,heigh=12)
  
}
EndothelialCells<-data.frame(Features=c("Pecam1","Cdh5"),CellType="EndothelialCells")
Artery<-data.frame(Features=c("Sox17","Dkk2","Efbn2","Sema3g","Fbln5","Hey1","Bmx","Mgp"),CellType="Artery")
Venous<-data.frame(Features=c("Ackr1","Selp","Sele","Icam1","Vwf","Vcam1","Nr2f2","Ephb4"),CellType="Venous")
Capillary<-data.frame(Features=c("Rgcc","Cd36","Kdr","Car2","Cd200","Cd300lg","Giphbp1","Aqp1"),CellType="Capillary")
TipCells<-data.frame(Features=c("Esm1","Dll4","Cxcr4","Col4a2","Col4a1","Sparc","Insr","Angpt2"),CellType="TipCells")
Lymphatics <- data.frame(Features=c("Prox1","Lyve1","Pdpn","Cd63","Flt4","Ccdc3","Mmrn2","Foxp2","Cldn11","Net1","Neo1"),CellType="Lymphatics")
Fenestration <- data.frame(Features=c("Plvap","Plpp3","Igfbp3"),CellType="Fenestration")
Angiogenesis <- data.frame(Features=c("Angpt1","Pdgfc","Vegfc","Vegfb","Lrp1","Tek","Kdr","Angpt2"),CellType="Angiogenesis")
Immature <- data.frame(Features=c("Pdlim1","Rplp0","Gapdh","Igfbp7","Rps2","Rpsa","Rpl12","Aplnr"),CellType="Immature")

Lista_Markers <- list(EndothelialCells=EndothelialCells,Artery=Artery,Venous=Venous,Capillary=Capillary,TipCells=TipCells,Lymphatics=Lymphatics,Fenestration=Fenestration,Angiogenesis=Angiogenesis,Immature=Immature)


source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_Annotations_Cluster.R")
source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_SeuratPipeline.R")

###############################################################################

data <- readRDS("0.2_SeuratPipeline/Seu.Obj.rds")
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

integration.var<-"ID"

dim.final <- Calculo.Dimensiones.PCA(data)
ElbowPlot(data, ndims = 30, reduction = "pca")
ggsave(paste(path.guardar,"Harmony_ElbowPlot_Log.png",sep="/"),width = 5,height = 5)

# Integration with Harmony
data.harmony <- RunHarmony(data,group.by.vars = integration.var,dims = 1:dim.final,plot_convergence=TRUE,reduction.save = "HarmonyLog",assay = "RNA")
data.harmony <- RunUMAP(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final)
data.harmony <- FindNeighbors(data.harmony,reduction = "HarmonyLog",dims = 1:dim.final,graph.name = "Harmony_Log")
resolutions <- c(0.1,0.3,0.5,0.7,0.9,1)
data.harmony <- FindClusters(data.harmony, resolution = resolutions,graph.name = "Harmony_Log")

# Save the data
print("Savind Data")
saveRDS(data.harmony,paste(path.guardar,"IntegratedData.rds",sep="/"))
print("Data Saved")

# 3)Plot the UMAPS and save them
col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.1)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.1.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.3)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.3", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.3.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.5)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.5", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.5.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.7)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.7", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.7.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.0.9)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.0.9", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.0.9.png",sep="/"),width=10,height=10)

col <-  getPalette(length(unique(data.harmony$Harmony_Log_res.1)))
DimPlot(data.harmony, group.by = "Harmony_Log_res.1", pt.size = 1,label = T,cols=col,raster=FALSE)&NoAxes()
ggsave(filename=paste(path.guardar,"LogNormalized_Harmony_res.1.png",sep="/"),width=10,height=10)

col <-  c("#99C5E3","#8CCE7D","#FFC685","#E2B5D5","#AAAEB0","#FA8D76","#F4D166")
DimPlot(data.harmony,reduction = "umap",group.by = "ID",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"ID.png",sep="/"),width = 10,height = 10)

col <-  c("#95BDED","#F5A040","#C8C2F2")
DimPlot(data.harmony,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE,split.by = "Phenotype")&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype_spli.png",sep="/"),width = 15,height = 7)

col <-  c("#95BDED","#F5A040","#C8C2F2")
DimPlot(data.harmony,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype.png",sep="/"),width = 7,height = 7)

col <-  c("#99C5E3","#8CCE7D","#FFC685","#E2B5D5","#AAAEB0","#FA8D76","#F4D166","#CFCAF5")
matchSCore2::summary_barplot(class.fac = data.harmony$Harmony_Log_res.0.3,obs.fac =data.harmony$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res03.png",sep="/"),width = 5,height = 7)

matchSCore2::summary_barplot(class.fac = data.harmony$Harmony_Log_res.0.5,obs.fac =data.harmony$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res05.png",sep="/"),width = 5,height = 7)

col <-  getPalette(20)
path.guardar <- "0.2_SeuratPipeline/Integration/Annotations"
dir.create(path.guardar)
path.guardar <- "0.2_SeuratPipeline/Integration/Annotations/Res0.3"
dir.create(path.guardar)
data.harmony <- SetIdent(data.harmony,value="Harmony_Log_res.0.3")
Annotation_Function(data.harmony,Lista_Markers,path.guardar)
markers <- FindAllMarkers(data.harmony,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res03.xlsx",sep="/"))

path.guardar <- "0.2_SeuratPipeline/Integration/Annotations/Res0.5"
dir.create(path.guardar)
data.harmony <- SetIdent(data.harmony,value="Harmony_Log_res.0.5")
Annotation_Function(data.harmony,Lista_Markers,path.guardar)
markers <- FindAllMarkers(data.harmony,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res05.xlsx",sep="/"))

path.guardar <- "0.2_SeuratPipeline/Integration/Annotations/Res0.7"
dir.create(path.guardar)
data.harmony <- SetIdent(data.harmony,value="Harmony_Log_res.0.7")
Annotation_Function(data.harmony,Lista_Markers,path.guardar)
markers <- FindAllMarkers(data.harmony,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res07.xlsx",sep="/"))





