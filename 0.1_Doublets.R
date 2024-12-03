# SCRIPT: Doublet Finder
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 6-11-2023

####################################################################################################################

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

# Load the script with the functions 
source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_DoubletDetection.R")
source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_SeuratPipeline.R")

#####################################################################################################################

# Load the data 
data<-readRDS("0.2_SeuratPipeline/Data.Filtered.rds")

path.guardar <-"0.2_SeuratPipeline/Doublets"
dir.create(path.guardar)

data.patient<-SplitObject(data,split.by="ID") #Object by patient

Patient.Data.Process <- foreach::foreach(i=1:length(data.patient),.final=function(x)setNames(x,names(data.patient)))%do%{
  print(paste("SeuratPipeline in",names(data.patient)[i],sep=" "))
  patient.data<-data.patient[[i]]
  patient.data<-SeuratPipeline(patient.data)}

# Doublet Estimation

Doublets.Estimation <- foreach::foreach(i=1:length(Patient.Data.Process),.final=function(x)setNames(x,names(Patient.Data.Process)))%do%{
  print(paste("Estimating Doublets",names(Patient.Data.Process)[i],sep=" "))
  data.patient<-Patient.Data.Process[[i]]
  data.patient<-DoubletDetection_DF(data.patient)
}

saveRDS(Doublets.Estimation,paste(path.guardar,"DoubletsRes.rds",sep="/")) 
