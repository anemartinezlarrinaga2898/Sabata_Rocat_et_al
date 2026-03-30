# SCRIPT: Doublet Finder
# AUTHOR: ANE MARTINEZ LARRINAGA
# DATE: 19-12-2023

# In this script we are going to use Doublet Finder to perfom the analysis. Doublet Finder is based in 
# generating doublets based on the average expression of two random cells. Afterwards what it does it
# runs the seurat preprocessing pipeline and in this case it computes a PCA reduction in order to generate
# a dimension where the artifical and the unreal data can be plotted and see which ones are more similar between each
# other. Then it introduces the metric that estimates that the real and unreal point should co.localize together. 

#####################################################################################################################

# Setting working parameters
source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")
directory<-setwd(Lista_Paths_Main$MainPath)

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(DoubletFinder)

# Define varibales to divided the object 
patient.columns<-"ID"

# Load the script with the functions 

source(paste(Lista_Paths_Main$Path_Utils,"0.2_Util_DoubletDetection.R",sep="/"))
source(paste(Lista_Paths_Main$Path_Utils,"0.3_Util_SeuratPipeline.R",sep="/"))

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.1_SeuratPipeline/Doublets",sep="/")
dir.create(path.guardar,recursive=TRUE)
#####################################################################################################################

# | Load the data 

path_obj<-paste(path.guardar_original,"0.1_SeuratPipeline",sep="/")

data <- readRDS(paste(path_obj,"Data.Filtered.rds",sep="/"))

# | Split the data bypatients
print("Splitting Seurat Objects by Patients")
data.patient<-SplitObject(data,split.by=patient.columns) #Object by patient
  
# Seurat Pipeline for each of the patients
print("Performing Seurat Pipeline in each patient")
  
Patient.Data.Process <- foreach::foreach(i=1:length(data.patient),.final=function(x)setNames(x,names(data.patient)))%do%{
    print(paste("SeuratPipeline in",names(data.patient)[i],sep=" "))
    patient.data<-data.patient[[i]]
    patient.data<-SeuratPipeline(patient.data)}
    
# | Doublet Estimation
    
Doublets.Estimation <- foreach::foreach(i=1:length(Patient.Data.Process),.final=function(x)setNames(x,names(Patient.Data.Process)))%do%{
    print(paste("Estimating Doublets",names(Patient.Data.Process)[i],sep=" "))
    data.patient<-Patient.Data.Process[[i]]
    data.patient<-DoubletDetection_DF(data.patient)}
  
saveRDS(Doublets.Estimation,paste(path.guardar,"DoubletFinder_ByPatient_DataFrame.rds",sep="/")) 

# | Add the metadata info 
print("names of Doublets.Estimation")
names(Doublets.Estimation)
print("Adding metadata")
data <- readRDS(paste(path_obj,"Data.Filtered.rds",sep="/"))
data.patient<-SplitObject(data,split.by=patient.columns) #Object by patient
  
Lista.Patient.Doublets<-vector(mode="list",length=length(data.patient))
names(Lista.Patient.Doublets)<-names(data.patient)

names(data.patient)

for(i in seq_along(data.patient)){
    patient<-names(data.patient)[i]
    print(patient)
    d.p<-data.patient[[i]] # Seurat Object of each one patient

    idx.results.db<-which(names(Doublets.Estimation)==patient)
    db.res <- Doublets.Estimation[[idx.results.db]] # Results of the doublet for the corresponding patient
    d.p@meta.data <- cbind(d.p@meta.data,db.res)
    Lista.Patient.Doublets[[i]]<-d.p
}

Seurat.Object.Total <-merge(x =  Lista.Patient.Doublets[[1]],y= Lista.Patient.Doublets[2:length( Lista.Patient.Doublets)])
saveRDS(Seurat.Object.Total,paste(path.guardar,"Doublets_Seurat_Total.rds",sep="/")) 
