# SCRIPT: Pasar las signatures de Ana por el objeto completo
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 29.07.2024

################################################################################
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
library(scGate)
library(org.Mm.eg.db)
library(clusterProfiler)
library(scGate)

getPalette <-  colorRampPalette(brewer.pal(8, "Set1"))
col <-  getPalette(20)
path.guardar <- "0.2_SeuratPipeline/Annotations"
dir.create(path.guardar)

################################################################################

data <- readRDS("0.2_SeuratPipeline/Seu.Obj_Remove.rds")
data$Phenotype <- factor(data$Phenotype,levels = c("OB_et","OB_40HT"),labels = c("OB_et","OB_40HT"))
data <- SetIdent(data,value = "RNA_snn_res.0.3")
# Annotate everything as signature .............................................

Signatures <- read_excel("0.0_Summaries/Signatures.xlsx")
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
  Vln<-VlnPlot(data, features =sig_names ,sort = T,cols=col)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
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
  
  Vln<-VlnPlot(data, features =sig_names ,sort = T,cols=col)& theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  file_name<-paste("Vln_",sig_names,".png",sep="")
  ggsave(plot=Vln,filename=paste(path.guardar,file_name,sep="/"),width=5,heigh=5)
}

# ..............................................................................

# Individual genes .............................................................

Lista_Signature_Complete <- list()

for(i in seq_along(Lista_Signature)){
  sig <- Lista_Signature[[i]]
  sig_name <- names(Lista_Signature)[i]
  print(sig_name)
  path.guardar <- "0.2_SeuratPipeline/Annotations"
  path.guardar <- paste(path.guardar,sig_name,sep="/")
  j <- 1
  for(j in seq_along(sig)){
    x <- sig[j]
    #print(x)
    idx <- which(rownames(data@assays$RNA@counts)==x)
    
    if(length(idx)>0){
      Fp<-FeaturePlot(data, features = x, min.cutoff = "q9", order = T) &NoAxes()
      file_name <- paste(x,".png",sep="")
      ggsave(plot=Fp,filename=paste(path.guardar,file_name,sep="/"),width=5,heigh=5)
    }else{
      next
    }
  }
}


# Division estimation 
path.guardar <- "0.2_SeuratPipeline/Annotations"
s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <-  stringr::str_to_title(cc.genes$g2m.genes)

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

col <-  c("#95BDED","#F5A040","#C8C2F2")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE,split.by = "Phase")&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype_SplitByPhase.png",sep="/"),width = 15,height = 7)

col <-  c("#95BDED","#F5A040","#C8C2F2")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype.png",sep="/"),width = 7,height = 7)

DimPlot(data,reduction = "umap",group.by = "Phase",label = TRUE, label.size = 5,cols =c("#95BDED","#F5A040","#C8C2F2"),pt.size = 0.5,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phase.png",sep="/"),width = 7,height = 7)

FeaturePlot(data, features = "S.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score.png",sep="/"),width = 7,height = 7)
FeaturePlot(data, features = "G2M.Score", min.cutoff = "q9", order = T) &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score.png",sep="/"),width = 7,height = 7)

FeaturePlot(data, features = "S.Score", min.cutoff = "q9", order = T,split.by = "Phenotype") &NoAxes()
ggsave(filename = paste(path.guardar,"S.Score_SplitByPhase.png",sep="/"),width = 10,height = 7)
FeaturePlot(data, features = "G2M.Score", min.cutoff = "q9", order = T,split.by = "Phenotype") &NoAxes()
ggsave(filename = paste(path.guardar,"G2M.Score_SplitByPhase.png",sep="/"),width = 10,height = 7)


