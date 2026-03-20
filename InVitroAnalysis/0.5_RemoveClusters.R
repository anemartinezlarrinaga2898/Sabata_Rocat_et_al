# SCRIPT: emove some of the clusters as they are noise or not clear what they are 
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 10.07.2024

###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library(Seurat)
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

data <- readRDS("0.2_SeuratPipeline/Seu.Obj.rds")
data <- SetIdent(data,value="RNA_snn_res.0.3")

data_cluster <- subset(data,idents = c("6","7"),invert=T)
gene.info.distribution <- summary(Matrix::colSums(data_cluster@assays$RNA@counts[,]>0))
hvg.number <- round(gene.info.distribution[4]+100)
data_cluster <- FindVariableFeatures(object = data_cluster,selection.method = "vst", nfeatures = hvg.number)
data_cluster <- ScaleData(object = data_cluster)

data_cluster <- RunPCA(object = data_cluster)

# Estimate the PCA Dimmensions. 
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
dim.final <- min(co1, co2)

data_cluster <- FindNeighbors(object = data_cluster, dims = 1:dim.final)
data_cluster <- RunUMAP(object = data_cluster, dims = 1:dim.final)
resolutions <- c(0.1,0.3,0.5)
data_cluster <- FindClusters(data_cluster, resolution = resolutions)

col <-  getPalette(length(unique(data_cluster$RNA_snn_res.0.5)))
DimPlot(data_cluster,reduction = "umap",group.by = "RNA_snn_res.0.1",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.1_Remove.png",sep="/"),width = 10,height = 10)

DimPlot(data_cluster,reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.3_Remove.png",sep="/"),width = 10,height = 10)

DimPlot(data_cluster,reduction = "umap",group.by = "RNA_snn_res.0.5",label = TRUE, label.size = 5,cols =col,pt.size = 1,raster=FALSE)&NoAxes()&NoLegend()
ggsave(filename = paste(path.guardar,"RNA_snn_res.0.5_Remove.png",sep="/"),width = 10,height = 10)

col <-  c("#99C5E3","#8CCE7D","#FFC685","#E2B5D5","#AAAEB0","#FA8D76","#F4D166")
DimPlot(data_cluster,reduction = "umap",group.by = "ID",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE)&NoAxes()
ggsave(filename = paste(path.guardar,"ID_Remove.png",sep="/"),width = 10,height = 10)

col <-  c("#99C5E3","#8CCE7D","#FFC685")
DimPlot(data,reduction = "umap",group.by = "Phenotype",label = FALSE, label.size = 5,cols =alpha(col,0.66),pt.size = 1,raster=FALSE,split.by = "Phenotype")&NoAxes()
ggsave(filename = paste(path.guardar,"Phenotype_spli_Remove.png",sep="/"),width = 15,height = 7)


col <-  getPalette(10)
matchSCore2::summary_barplot(class.fac = data_cluster$RNA_snn_res.0.1,obs.fac =data_cluster$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res01_Remove.png",sep="/"),width = 5,height = 7)

matchSCore2::summary_barplot(class.fac = data_cluster$RNA_snn_res.0.3,obs.fac =data_cluster$Phenotype)+scale_fill_manual(values = col)
ggsave(filename = paste(path.guardar,"BarPlot_Pheno_Res03_Remove.png",sep="/"),width = 5,height = 7)

Fp <- FeaturePlot(data_cluster, features = c("Pecam1","Cdh5","Vwf"), min.cutoff = "q9", order = T,pt.size = 1,cols=c("Grey","Red"))&NoAxes()
ggsave(plot = Fp,filename = paste(path.guardar,"FeaturePlots_MarkersEndo_Remove.png",sep="/"),width = 10,height = 10)

Fp <- FeaturePlot(data_cluster, features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Capillary_Remove.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_cluster, features = c("Sox17","Hey1","Sema3g","Clu"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Artery_Remove.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_cluster, features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Venous_Remove.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_cluster, features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Tip_Remove.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_cluster, features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Division_Remove.png",sep="/"),plot=Fp,width=10,height=10)
Fp <- FeaturePlot(data_cluster, features = c("Lyve1","Prox1","Pdpln"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Lymphatis_Remove.png",sep="/"),plot=Fp,width=10,height=5)
Fp <- FeaturePlot(data_cluster, features = c("Lyve1","Prox1","Hey1","Hey2","Unc5b","Flt4","Nrp2","Nr2f2","Ephb4","Cxcr4"), min.cutoff = "q9", order = T,raster = FALSE,cols=c("Grey","Red"))&NoAxes()
ggsave(filename=paste(path.guardar,"Fp_Interes_Remove.png",sep="/"),plot=Fp,width=15,height=10)

col <-  getPalette(10)
data_cluster <- SetIdent(data_cluster,value="RNA_snn_res.0.3")
Vln <- VlnPlot(object = data_cluster,features = c("Kdr","Rgcc","Cd200","Cd300lg","Cd36","Sgk1"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Capillary_Remove.png",sep="/"),plot=Vln,width=10,height=10)
Vln <- VlnPlot(object = data_cluster,features = c("Sox17","Hey1","Sema3g","Clu"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Artery_Remove.png",sep="/"),plot=Vln,width=10,height=10)
Vln <- VlnPlot(object = data_cluster,features = c("Nr2f2","Vcam1","Vwf","Vcam1","Icam1"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Venous_Remove.png",sep="/"),plot=Vln,width=10,height=10)
Vln <- VlnPlot(object = data_cluster,features = c("Esm1","Cxcr4","Dll4","Col4a1","Col4a2"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Tip_Remove.png",sep="/"),plot=Vln,width=10,height=10)
Vln <- VlnPlot(object = data_cluster,features = c("Mki67","Cdk1","Cdk2","Cdk4","Cdk6"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Division_Remove.png",sep="/"),plot=Vln,width=10,height=10)
Vln <- VlnPlot(object = data_cluster,features = c("Lyve1","Prox1","Pdpln"),cols = col,pt.size = 0.1,sort = T)
ggsave(filename=paste(path.guardar,"Vln_Lymphatics_Remove.png",sep="/"),plot=Vln,width=10,height=10)



col <-  getPalette(length(unique(data_cluster$RNA_snn_res.0.5)))
data_cluster <- SetIdent(data_cluster,value="RNA_snn_res.0.1")
saveRDS(data_cluster,paste(path.guardar,"Seu.Obj_Remove.rds",sep="/"))

data_cluster <- SetIdent(data_cluster,value="RNA_snn_res.0.1")
markers <- FindAllMarkers(data_cluster,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res01_Remove.xlsx",sep="/"))
markers_cluster <- split(markers,markers$cluster)
openxlsx::write.xlsx(markers_cluster,paste(path.guardar,"Markers_ByCluster_01_Remove.xlsx",sep="/"))


data_cluster <- SetIdent(data_cluster,value="RNA_snn_res.0.3")
markers <- FindAllMarkers(data_cluster,only.pos = T,min.pct = 0.25)
writexl::write_xlsx(markers,paste(path.guardar,"TotalMarkers_Res03_Remove.xlsx",sep="/"))
markers_cluster <- split(markers,markers$cluster)
openxlsx::write.xlsx(markers_cluster,paste(path.guardar,"Markers_ByCluster_03_Remove.xlsx",sep="/"))



