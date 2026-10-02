################################################################################
# SCRIPT: DoubletFinder estimation for MUT samples
# AUTHOR: Ane Martinez Larrinaga
# DATE: 19-12-2023
#
# DESCRIPTION:
# Runs the Seurat preprocessing workflow independently for each MUT sample, estimates doublets using DoubletFinder, stores the predictions in sample metadata, and merges the samples. Predicted doublets are retained and are not removed at this stage.
#
# INPUT:
#   results/0.5_MUT_Analysis/Obj/Data.Filtered.rds
#
# OUTPUT:
#   results/0.5_MUT_Analysis/Total/Doublets/Doublets_Seurat_Total.rds
#
# NOTE:
#   Analytical thresholds, clustering resolutions and biological selections
#   are preserved from the original analysis unless explicitly documented in
#   REVIEW_NOTES.md.
################################################################################

# ------------------------------------------------------------------------------
# Project-relative paths
# ------------------------------------------------------------------------------

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = "."),
  mustWork = FALSE
)
results_root <- file.path(project_root, "results")
utils_dir <- file.path(project_root, "utils")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(DoubletFinder)

# Metadata field used to split the object 
patient.columns <- "ID"

# Load project utility functions

source(file.path(utils_dir, "0.2_Util_DoubletDetection.R"))
source(file.path(utils_dir, "0.3_Util_SeuratPipeline.R"))

# Define output paths
path.guardar_original <- results_root
path.guardar<-paste(path.guardar_original,"0.5_MUT_Analysis/Total/Doublets",sep="/")
dir.create(path.guardar,recursive=TRUE)
################################################################################

# Load input data

path_obj<-paste(path.guardar_original,"0.5_MUT_Analysis/Obj",sep="/")

data <- readRDS(paste(path_obj,"Data.Filtered.rds",sep="/"))

# | Split the data by samples
print("Splitting Seurat objects by sample")
data.patient<-SplitObject(data,split.by=patient.columns) # Object by sample
  
# Run the Seurat pipeline independently for each sample
print("Running the Seurat pipeline for each sample")
  
Patient.Data.Process <- foreach::foreach(i=1:length(data.patient),.final=function(x)setNames(x,names(data.patient)))%do%{
    print(paste("SeuratPipeline in",names(data.patient)[i],sep=" "))
    patient.data<-data.patient[[i]]
    patient.data<-SeuratPipeline(patient.data)}
    
# Estimate doublets
    
Doublets.Estimation <- foreach::foreach(i=1:length(Patient.Data.Process),.final=function(x)setNames(x,names(Patient.Data.Process)))%do%{
    print(paste("Estimating Doublets",names(Patient.Data.Process)[i],sep=" "))
    data.patient<-Patient.Data.Process[[i]]
    data.patient<-DoubletDetection_DF(data.patient)}
  
saveRDS(Doublets.Estimation,paste(path.guardar,"DoubletFinder_ByPatient_DataFrame.rds",sep="/")) 

# Add DoubletFinder results to sample metadata
print("names of Doublets.Estimation")
names(Doublets.Estimation)
print("Adding metadata")
data <- readRDS(paste(path_obj,"Data.Filtered.rds",sep="/"))
data.patient<-SplitObject(data,split.by=patient.columns) # Object by sample
  
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
