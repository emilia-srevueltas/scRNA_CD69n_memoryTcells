#After file '1_Subsetting_sixpopulations', CD4 and CD8 cells were analyzed separately

library(SeuratObject)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(eulerr)
library(tidyverse)
library(patchwork)
library("BiocParallel")

setwd("/Users/EmiliaSR/Desktop/scRNA/Emilia")
getwd()

file_name = "CD4_29012024.RData"

load(file_name)
Idents(CD4)

#Author: Pawel Durek, adapted by Emilia Schneider Revueltas
#Cluster WITH integration: 

SO_all <- CD4

splitingGroup = "experiment"

DefaultAssay(SO_all) <- "RNA"


######################################################

#

# Integration

#

######################################################

SubsetSO.list <- SplitObject(SO_all, split.by = "experiment")

#rm(SO_all)


sublist.analysis.wrapper <- function(i){
  
  print(paste0("analyse experiment ", i))
  
  x <- NormalizeData(SubsetSO.list[[i]], verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 50, verbose = FALSE)
  
  return(x)
  
}

length(SubsetSO.list)
SubsetSO.list <- bplapply(1:length(SubsetSO.list), 
                          sublist.analysis.wrapper, 
                          BPPARAM=MulticoreParam(2))


SO.anchors    <- FindIntegrationAnchors(object.list = SubsetSO.list, 
                                        reduction = "rpca", dims = 1:50)
SO.integrated <- IntegrateData(anchorset = SO.anchors, dims = 1:50)


SO.integrated <- ScaleData(SO.integrated, verbose = FALSE)
SO.integrated <- RunPCA(SO.integrated, npcs = 50, verbose = FALSE)

SO.integrated <- RunUMAP(SO.integrated, reduction = "pca", dims = 1:50)
SO.integrated <- FindNeighbors(SO.integrated, dims = 1:50, reduction = "pca")



#integrated_snn_res.1

SO.integrated <- FindClusters(SO.integrated, resolution = 1, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.9, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.8, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.7, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.6, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.5, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.4, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.3, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.2, reduction = "pca")
SO.integrated <- FindClusters(SO.integrated, resolution = 0.1, reduction = "pca")


DimPlot(SO.integrated)

### Organize MetaData ###
SO.integrated@meta.data$tissue = SO.integrated@meta.data$experiment
SO.integrated@meta.data$tissue = as.factor(SO.integrated@meta.data$tissue)
levels(SO.integrated@meta.data$tissue) = c(rep("Blood",3),rep("BM",3))
levels(SO.integrated@meta.data$tissue)

SO.integrated@meta.data$donor = SO.integrated@meta.data$experiment
SO.integrated@meta.data$donor = as.factor(SO.integrated@meta.data$donor)
levels(SO.integrated@meta.data$donor)
levels(SO.integrated@meta.data$donor) = c("D1", "D2", "D3", "D1", "D2", "D3")

SO.integrated@meta.data$celltype = as.factor(SO.integrated@meta.data$celltype)
levels(SO.integrated@meta.data$celltype)

setwd("/Users/EmiliaSR/Desktop/scRNA/R Stuff/CD4")
save(SO.integrated, file = "CD4_integrated_29012024.RData")

Idents(SO.integrated)<- "celltype"
Idents(SO.integrated)<- "integrated_snn_res.0.2"

DimPlot(SO.integrated, split.by = "donor")
DimPlot(SO.integrated, split.by = "celltype")

CD4_celltypeFAM <- FindAllMarkers(SO.integrated, only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

write.csv(CD4_celltypeFAM, 
          file ="/Users/EmiliaSR/Desktop/scRNA/R stuff/CD4/CD4_FAM_celltype.csv")

####Quality control####

DefaultAssay(SO.integrated) <- "RNA"
SO.integrated[["percent.mito"]] <- PercentageFeatureSet(object = SO.integrated, pattern = "^MT-")

VlnPlot(SO.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
