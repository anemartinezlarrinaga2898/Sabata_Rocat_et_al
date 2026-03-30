# SCRIPT: Signature Estimation Using Add Module Score
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 25.03.2025

#######################################################################################################################

source("/ijc/USERS/amartinezl/MGRAUPERA_18/Paths.R")

directory<-setwd(Lista_Paths_Main$MainPath)


library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(foreach)
library(harmony)
library(Rcpp)
library(readxl)
library(DropletUtils)
library(DelayedArray)

# Obtener path de anotaciones endoteliales
path.guardar_original <- Lista_Paths_Main$Path_Guardar_Retinas
path.guardar<-paste(path.guardar_original,"0.1_SeuratPipeline",sep="/")
dir.create(path.guardar,recursive=TRUE)
path_obj<-"/ijc/LABS/GRAUPERA/DATA/ANE_MARTINEZ_LARRINAGA/MGRAUPERA_18_ARI/Nov_2025"

# | QC definition for future QC control

pattern.mito<-"^mt-"
pattern.ribo<-"^rps"

# Function 

Read_CellBender_h5_Mat <- function(file_name,use.names = TRUE,unique.features = TRUE) {
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }

  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }

  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")

  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]


  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )

  if (unique.features) {
    features <- make.unique(names = features)
  }

  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")

  infile$close_all()

  return(sparse.mat)
}
###################################################################################################

# | List Files 

dirs <- list.dirs(path_obj, recursive = FALSE, full.names = FALSE)
dirs_validas <- dirs[grepl("_output$", dirs) & !grepl("_old$", dirs)]

sample_names<-dirs_validas
samples_names<-stringr::str_replace_all(sample_names,"_output","")

Lista_Seurat_Object<-vector(mode="list",length=length(samples_names))
names(Lista_Seurat_Object)<-samples_names

# | Load the Seurat Object
for(i in seq_along(dirs_validas)){
    sn<-dirs_validas[i]
    sn_name<-stringr::str_replace_all(sn,"_output","")
    print(sn)

    pheno<-sample_pheno[which(names(sample_pheno) == sn_name)]
    names(pheno)<-NULL
    path_folder<-paste(path_obj,sn,"outs/cellbender",sep="/")
    file_name<-paste(sn_name,"_cellbender_filtered.h5",sep="")
    counts<-Read_CellBender_h5_Mat(file.path(path_folder, file_name))

    SeuratObject <- CreateSeuratObject(counts = counts,project = "Retinas_eGFP")
    SeuratObject$ID<-sn_name
    SeuratObject$Phenotype<-pheno

    SeuratObject$log10GenesPerUMI <- log10(SeuratObject$nFeature_RNA) / log10(SeuratObject$nCount_RNA)
    SeuratObject$percent.mt <- PercentageFeatureSet(SeuratObject, pattern = pattern.mito)
    SeuratObject$percent.mt.div100 <- SeuratObject$percent.mt/100
    SeuratObject$percent.rb <- PercentageFeatureSet(SeuratObject, pattern = pattern.ribo)

    Lista_Seurat_Object[[i]]<-SeuratObject
}

Seurat.Object.Total <-merge(x = Lista_Seurat_Object[[1]],y=Lista_Seurat_Object[2:length(Lista_Seurat_Object)])
Seurat.Object.Total <- SetIdent(Seurat.Object.Total,value="ID")
saveRDS(Seurat.Object.Total,paste(path.guardar,"Seurat.Object.Total.rds",sep="/"))

# | QC Plotting

# - ViolinPlot total 
VlnPlot(Seurat.Object.Total, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC.png",sep="/"),width = 15,height = 8)

# - Complexity 
Complexity <- Seurat.Object.Total@meta.data %>% ggplot(aes(x=log10GenesPerUMI, color = Phenotype, fill=Phenotype)) +
    geom_density(alpha = 0.2) +
    theme_classic()
ggsave(filename = paste(path.guardar,"Complexity.png",sep="/"),plot = Complexity,width = 10,height = 10)

# - CellCounts
print("Plotting CellCounts")
CellCounts <- Seurat.Object.Total@meta.data %>% ggplot(aes(x=ID, fill=ID)) +
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
ggsave(filename = paste(path.guardar,"CellCounts.png",sep="/"),plot = CellCounts,width = 35,height = 10)

# - UMICount
print("Plotting CellCounts")
UMICount <-Seurat.Object.Total@meta.data %>%
    ggplot(aes(color=ID, x=nCount_RNA, fill= ID)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
ggsave(filename = paste(path.guardar,"UMICount.png",sep="/"),plot = UMICount,width = 35,height = 10)

# - GenesPerCell
print("Plotting GenesPerCell")
GenesPerCell <-Seurat.Object.Total@meta.data %>% ggplot(aes(color=ID, x=nFeature_RNA, fill= ID)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_x_log10() +
    geom_vline(xintercept = 300)
ggsave(filename = paste(path.guardar,"GenesPerCell.png",sep="/"),plot = GenesPerCell,width = 35,height = 10)


# | Filtering samples
data.filtered <- subset(Seurat.Object.Total,subset=(nFeature_RNA>=250) & (percent.mt<20) & (log10GenesPerUMI>0.8))
VlnPlot(data.filtered, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"AfterQC.png",sep="/"),width = 15,height = 8)
saveRDS(data.filtered,paste(path.guardar,"Data.Filtered.rds",sep="/"))
