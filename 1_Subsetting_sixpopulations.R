#Upstream Analysis of the scRNA of ex-vivo resting human memory T cells from paired blood and bone marrow: 
#demultiplexing, mapping, detection of intact cells and quantification of gene expression. GEX-libraries were 
#merged without normalization and re analyze without TCR related genes. 
# Upstream Analysis was performed by Pawel Durek and Frederik Heinrich. 


library(SeuratObject)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(eulerr)
library(tidyverse)
library(patchwork)
library("BiocParallel")

setwd("/Users/EmiliaSR/Desktop/scRNA/R Stuff")


file_name = "Emila_noXCR.integrated_SNN_rpca_50.RData"
load(file_name)

DefaultAssay(SO.integrated) <- "RNA"

#Manual gating on scattered log-normalized CITE-seq counts was used to define 6 populations:
#'Blood_CD4p_CD69n','BM_CD4p_CD69n', 'BM_CD4p_CD69p', 'Blood_CD8p_CD69n','BM_CD8p_CD69n', 'BM_CD8p_CD69p',
#where  CD4p corresponds to CITE-seq counts of CD4+ cells, CD8p to CD8+ cells and CD69p to surface CD69+ and CD69n to CD69-.
#This process corresponds to Supplementary Figure 1c of the Manuscript. 

####Add annotations from Scattered gating antibody based####
setwd("/Users/EmiliaSR/Desktop/scRNA/Emilia")

annotation_file <- "CD69_CD8_4_split.cloupe.csv"
annot2 <- read.table(annotation_file, sep = ",", row.names = c(1), header = TRUE)
annot2[,2] <- NULL
annot2[,1] <- NULL

SO.integrated <- AddMetaData(SO.integrated, metadata = annot2)
print(colnames(SO.integrated@meta.data))
print(unique(SO.integrated@meta.data$Groups_Lib))

SO.integrated@meta.data$Groups_Lib = as.factor(SO.integrated@meta.data$Groups_Lib)
levels(SO.integrated@meta.data$Groups_Lib)

SO.integrated@meta.data$celltype = SO.integrated@meta.data$Groups_Lib
levels(SO.integrated@meta.data$celltype)

####Subset desired population as new seurat objet####

CD4 <- subset(x = SO.integrated, subset = (celltype == 'Blood_CD4p_CD69n'
                                          | celltype == 'BM_CD4p_CD69n'
                                          | celltype == 'BM_CD4p_CD69p' ))
unique(CD4@meta.data$celltype)
Idents(CD4) <- "celltype"
View(CD4)
Idents(CD4)
save(CD4, file = "CD4_29012024.RData")


CD8 <- subset(x = SO.integrated, subset = (celltype == 'Blood_CD8p_CD69n'
                                           | celltype == 'BM_CD8p_CD69n'
                                           | celltype == 'BM_CD8p_CD69p' ))
unique(CD8@meta.data$celltype)
Idents(CD8) <- "celltype"
View(CD8)
Idents(CD8)
save(CD8, file = "CD8_29012024.RData")

