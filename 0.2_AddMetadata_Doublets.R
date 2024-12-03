# SCRIPT: Doublet Finder
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 10-10-2023

#####################################################################################################################

# Setting working parameters

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

getPalette <-  colorRampPalette(brewer.pal(9, "Blues"))
colors <-  sample(getPalette(11))

library(patchwork)
library(foreach)
library(DoubletFinder)
library(yaml)

# Define working paths
directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

#####################################################################################################################

# 0) Load the DoubletFinder Results

doublet.finder.results<-readRDS("0.2_SeuratPipeline/Doublets/DoubletsRes.rds") 
# Load the Seurat Object
data<-readRDS("0.2_SeuratPipeline/Data.Filtered.rds")
data.patient<-SplitObject(data,split.by="ID")  #Object by patient

Lista.Patient.Doublets<-vector(mode="list",length=length(data.patient))
names(Lista.Patient.Doublets)<-names(data.patient)

for(i in seq_along(data.patient)){
  patient<-names(data.patient)[i]
  d.p<-data.patient[[i]] # Seurat Object of each one patient
  
  idx.results.db<-which(names(doublet.finder.results)==patient)
  db.res <- doublet.finder.results[[idx.results.db]] # Results of the doublet for the corresponding patient
  d.p@meta.data <- cbind(d.p@meta.data,db.res)
  Lista.Patient.Doublets[[i]]<-d.p
}

Seurat.Object.Total <-merge(x =  Lista.Patient.Doublets[[1]],y= Lista.Patient.Doublets[2:length( Lista.Patient.Doublets)])

saveRDS(Seurat.Object.Total,"0.2_SeuratPipeline/Data_Doublets.rds") 