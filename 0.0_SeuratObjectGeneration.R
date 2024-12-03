# SCRIPT: Generate Seurat Object by sample
# AUTOR: ANE MARTINEZ LARRINAGA
# FECHA: 23-10-2023

# has creado el objeto en el cluster porque aqui te lo crea con la version V5
###############################################################################

directory <- setwd("/Users/anemartinezlarrinaga/Library/CloudStorage/OneDrive-JosepCarrerasLeukaemiaResearchInstitute(IJC)/2_PhD/1_GRAUPERA_LAB/2_PROYECTOS/10_Ana_MGRAUPERA_10/")

library("Seurat", lib.loc = "/Library/Frameworks/R.framework/Versions/4.2/old-versions")
options(Seurat.object.assay.version = 'v4')
library(RColorBrewer)
library(tidyverse)
library(foreach)
library(Matrix)
library(readxl)

getPalette <-  colorRampPalette(brewer.pal(9, "Paired"))
col <-  getPalette(10)
options(Seurat.object.assay.version = 'v4')
path.guardar <- "0.2_SeuratPipeline"
dir.create(path.guardar)
###############################################################################

files.to.read <- list.files("0.1_CellRanger")

# 0) Read the matrix
Files <- foreach (i =1:length(files.to.read),.final =function(x)setNames(x,files.to.read))%do%{
  files <- files.to.read[i]
  print(files)
  files_path <- paste("0.1_CellRanger",files,sep="/")
  count_matrix <- Read10X(files_path)
}

# 1) Generate seurat object by mouse 

metadata <- readxl::read_excel("0.0_Summaries/SamplesInfo.xlsx")

Seu.QC <- vector(mode="list",length = length(Files))
names(Seu.QC) <- names(Files)

pattern.mito<-"^mt-"
pattern.ribo<-"^rps"

for(i in seq_along(Files)){
  cm <- Files[[i]]
  sample <- names(Files)[i]
  sample <- stringr::str_replace_all(sample,"_output","")
  print(sample)
  metadata.s <- metadata[which(metadata$Sample.Name==sample),]
  
  data <- CreateSeuratObject(cm,project = "mgraupera_10")
  data$ID <- metadata.s$Sample.Name
  data$Phenotype <- metadata.s$Phenotype
  
  data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA)
  data$percent.mt <- PercentageFeatureSet(data, pattern = pattern.mito)
  data$percent.mt.div100 <- data$percent.mt/100
  data$percent.rb <- PercentageFeatureSet(data, pattern = pattern.ribo)
  Seu.QC[[i]] <- data
}

data.total <-merge(x = Seu.QC[[1]],y=Seu.QC[2:length(Seu.QC)],add.cell.ids = names(Files))

data.total <- SetIdent(data.total,value="ID")
saveRDS(data.total,paste(path.guardar,"Data.Total.rds",sep="/"))
VlnPlot(data.total, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"QC.png",sep="/"),width = 15,height = 8)

data.filtered <- subset(data.total,subset=(nFeature_RNA>=250) & (percent.mt<20) & (log10GenesPerUMI>0.8))
VlnPlot(data.filtered, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, log = T,raster=FALSE)
ggsave(filename = paste(path.guardar,"AfterQC.png",sep="/"),width = 15,height = 8)

saveRDS(data.filtered,paste(path.guardar,"Data.Filtered.rds",sep="/"))


