# SCRIPT: Generate Seurat Object by sample
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 23-10-2023

###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library("Seurat", lib.loc = "/Library/Frameworks/R.framework/Versions/4.2/old-versions")
options(Seurat.object.assay.version = 'v4')
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)

getPalette <-  colorRampPalette(brewer.pal(8, "Set1"))
col <-  getPalette(10)
path.guardar <- "0.2_SeuratPipeline"

#source("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/3_UTILS/Util_Annotations_Cluster.R")
###############################################################################

data <- readRDS("0.2_SeuratPipeline/Data_Doublets.rds")

# Filtering genes in order to avaoid small fragments
n_cells <- Matrix::rowSums(data@assays$RNA@counts > 0)
kept_genes <- rownames(data)[n_cells > 5]
data <- subset(data, features = kept_genes)

data <- NormalizeData(object = data)
gene.info.distribution <- summary(Matrix::colSums(data@assays$RNA@counts[,]>0))
hvg.number <- round(gene.info.distribution[4]+100)
data <- FindVariableFeatures(object = data,selection.method = "vst", nfeatures = hvg.number)
data <- ScaleData(object = data)

data <- RunPCA(object = data)

# Estimate the PCA Dimmensions. 
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
dim.final <- min(co1, co2)


data <- FindNeighbors(object = data, dims = 1:dim.final)
data <- RunUMAP(object = data, dims = 1:dim.final)
resolutions <- c(0.1,0.3,0.5)
data <- FindClusters(data, resolution = resolutions)

col <-  getPalette(length(unique(data$RNA_snn_res.0.5)))
DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.1.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.3.png",sep="/"),width = 10,height = 10)

DimPlot(data,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.5.png",sep="/"),width = 10,height = 10)

col <-  c("#99C5E3","#8CCE7D","#FFC685","#E2B5D5","#AAAEB0","#FA8D76","#F4D166")
DimPlot(data,reduction = "umap",group.by = "ID",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"ID.png",sep="/"),width = 10,height = 10)

col <-  c("#99C5E3","#8CCE7D","#FFC685")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype.png",sep="/"),width = 10,height = 10)

col <-  c("#99C5E3","#8CCE7D","#FFC685")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE,split.by = "Phenotype")&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype_split.png",sep="/"),width = 15,height = 7)

col <-  getPalette(10)
matchSCore2::summary_barplot(class.fac = data$RNA_snn_res.0.1,obs.fac =data$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res01.png",sep="/"),width = 5,height = 7)

matchSCore2::summary_barplot(class.fac = data$RNA_snn_res.0.3,obs.fac =data$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res03.png",sep="/"),width = 5,height = 7)

FeaturePlot(data, features = c("Pecam1","Cdh5","Vwf"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(plot = x,filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

FeaturePlot(data, features = c("Pecam1","Cdh5","Vwf"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(plot = x,filename = paste(path.guardar,"FeaturePlots_MarkersEndo.png",sep="/"),width = 10,height = 10)

Fp <- FeaturePlot(data, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes.png",sep="/"),plot=Fp,width=15,height=10)


col <-  getPalette(length(unique(data$RNA_snn_res.0.5)))
data <- SetIdent(data,value="RNA_snn_res.0.1")
saveRDS(data,paste(path.guardar,"Seu.Obj.rds",sep="/"))

data <- SetIdent(data,value="RNA_snn_res.0.1")
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res01.xlsx",sep="/"))
markers_cluster <- split(markers,markers$cluster)
openxlsx::write.xlsx(markers_cluster,paste(path.guardar,"Markers_ByCluster.xlsx",sep="/"))


data <- SetIdent(data,value="RNA_snn_res.0.3")
markers <- FindAllMarkers(data,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res03.xlsx",sep="/"))
markers_cluster <- split(markers,markers$cluster)
openxlsx::write.xlsx(markers_cluster,paste(path.guardar,"Markers_ByCluster_03.xlsx",sep="/"))
