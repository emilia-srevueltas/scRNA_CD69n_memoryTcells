#From the integrated samples, here we removed cells lacking TCR information

library(Matrix)
library(usethis)
library(devtools)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(Rcpp)
library(harmony)

#install.packages("harmony")


setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8")

file_name = "CD8_integrated_29012024.RData"

load(file_name)

SO <- SO.integrated

####Pawel Freddy Colors ####

Disc_colors = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",
                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")

#### Choose IDENTS####
#Run only desired IDENT
Idents(SO) <- "integrated_snn_res.0.2"
DefaultAssay(SO) <- "RNA"

####Quality control####

DefaultAssay(SO) <- "RNA"
SO[["percent.mito"]] <- PercentageFeatureSet(object = SO, pattern = "^MT-")

VlnPlot(SO, pt.size= 0, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, cols = Disc_colors)


####Basic UMAP####

UMAP <- DimPlot(SO, pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP, id = "ident", size = 4, box = T, repel = F) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))


UMAP2 <- DimPlot(SO, split.by = "donor", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP2, id = "ident", size = 3) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP3 <- DimPlot(SO, split.by = "celltype", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP3, id = "ident", size = 3) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

U2_U3 <- (UMAP2 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))) / (UMAP3 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)))
U2_U3 

####Add TCRs####
setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/TCR_verified/vdj_t")

tcr_v <- read.csv("filtered_contig_annotations.csv")
tcr_v <- tcr_v[,c("barcode", "cdr3","cdr3_nt", "chain")]
tcr_v$TRA <- NA
tcr_v$TRB <- NA

tcr_df <- tcr_v %>%
  mutate(TRA = case_when(chain=="TRA" ~ as.character(cdr3))) %>%
  mutate(TRB = case_when(chain=="TRB" ~ as.character(cdr3)))

df_TRA <- tcr_df %>% filter(tcr_df$chain == "TRA")
df_TRB <- tcr_df %>% filter(tcr_df$chain == "TRB")

rownames(df_TRA) <- df_TRA[,1]
df_TRA[,1] <- NULL #I repeated until only names and TRA sequence stayed-> 4X
df_TRA[,1] <- NULL
df_TRA[,1] <- NULL
df_TRA[,1] <- NULL
df_TRA[,2] <- NULL

rownames(df_TRB) <- df_TRB[,1]
df_TRB[,1] <- NULL #I repeated until only names and TRB sequence stayed-> 5X
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL
df_TRB[,1] <- NULL

aa_cl <- merge(df_TRA, df_TRB, by = "row.names")
rownames(aa_cl) <- aa_cl[,1]
aa_cl[,1] <- NULL


aa_cl$TRA_TRB <- paste(aa_cl$TRA, aa_cl$TRB, sep = ":")
aa_cl[,1] <- NULL
aa_cl[,1] <- NULL


SO <- AddMetaData(object=SO, 
                  metadata=aa_cl, 
                  col.name = "aa_clones")

SO@meta.data$aa_clones = as.factor(SO@meta.data$aa_clones)
levels(SO@meta.data$aa_clones)


####PLOT CELLS WITHOUT TCRs ####

Plot_tmp = as.data.frame(SO@reductions$umap@cell.embeddings)

Plot_tmp$clonotype = SO@meta.data$aa_clones

Plot_tmp$col2 = SO@meta.data$integrated_snn_res.0.4
levels(SO@meta.data$integrated_snn_res.0.4)

Plot_tmp$col = ifelse(is.na(Plot_tmp$clonotype),"No clonotype",as.character(Plot_tmp$col2))
Plot_tmp$col <- factor(Plot_tmp$col, levels = c(levels(SO@meta.data$integrated_snn_res.0.4)))
use_these_colors = c(Disc_colors[1:length(levels(Plot_tmp$col2))],"#000000")
Plot = Plot_tmp[order(Plot_tmp$col),]

p <- ggplot(Plot, aes(UMAP_1, UMAP_2, color= col)) +
  geom_point(size=0.3) +
  labs(title= paste("UMAP"," highlighting ","cells without TCR",sep = ""),color="") +
  xlab("Umap 1") +
  ylab("Umap 2") +
  scale_color_manual(values = use_these_colors,na.value = "#000000") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title=element_text(size=20),
        legend.text=element_text(size=10,face="bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)) +
  guides(color = guide_legend(override.aes = list(size = 5)))
p

####Number of cells with TCR####

cell_with_tcrs <- which(!is.na(SO@meta.data$aa_clones))

cell_with_tcrs_barcodes <- row.names(SO@meta.data)[cell_with_tcrs]

Cell_Nr <- as.data.frame(summary(SO@meta.data$integrated_snn_res.0.2[cell_with_tcrs]))
Cell_Nr$totalCells <- summary(SO@meta.data$integrated_snn_res.0.2)
print(Cell_Nr)

write_csv(Cell_Nr, file ="/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8/CD8_CellNr_TCRs.csv")

Cell_Nr_donor <- as.data.frame(summary(SO@meta.data$donor[cell_with_tcrs]))
Cell_Nr_donor$totalCells <- summary(SO@meta.data$donor)
print(Cell_Nr_donor)

write_csv(Cell_Nr_donor, file ="/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8/CD8_CellNr_TCRs_donor.csv")

TCR_list_D1 <- as.factor((SO@meta.data$aa_clones)[SO@meta.data$donor == "D1"])
TCR_list_D2 <- as.factor((SO@meta.data$aa_clones)[SO@meta.data$donor == "D2"])
TCR_list_D3 <- as.factor((SO@meta.data$aa_clones)[SO@meta.data$donor == "D3"])

#setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4")
#write.table(cbind(Name=colnames(SO@meta.data), t(SO@meta.data)), 
            #file="Emilias_CD4_TCR_data.mwt.txt",append=FALSE, 
            #quote=FALSE,sep="\t",row.names = FALSE)
getwd()

####subset no TCRs####

summary(is.na(SO@meta.data$aa_clones))

# Convert aa_clones to character type
SO@meta.data$aa_clones <- as.character(SO@meta.data$aa_clones)

# Assign 'NoTCR' to missing values
SO@meta.data$aa_clones[is.na(SO@meta.data$aa_clones)] <- 'NoTCR'

# Convert aa_clones to character type
SO@meta.data$aa_clones <- as.factor(SO@meta.data$aa_clones)

# Check if it worked
summary(SO@meta.data$aa_clones[SO@meta.data$aa_clones == "NoTCR"])


SO_TCR <- subset(x = SO, subset = (aa_clones == "NoTCR"), invert = T)

summary(SO_TCR@meta.data$aa_clones[SO_TCR@meta.data$aa_clones == "NoTCR"])

summary(is.na(SO_TCR@meta.data$aa_clones))


####HARMONY INTEGRATION#### Author: Pawel Durek, Adapted: Emilia Schneider Revueltas

DefaultAssay(SO_TCR) 

SO_all <- SO_TCR

splitingGroup = "experiment"

DefaultAssay(SO_all) <- "RNA"

Idents(SO_all) <- "celltype" #Are Idents important?


SubsetSO.list <- SplitObject(SO_all, split.by = "experiment")

SO.harmony <- FindVariableFeatures(SO_all, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
SO.harmony <- ScaleData(SO.harmony, verbose = FALSE)
SO.harmony <- RunPCA(SO.harmony, npcs = 50, verbose = FALSE)
SO.harmony <- RunHarmony(SO.harmony, c(splitingGroup))

SO.harmony <- RunUMAP(SO.harmony, reduction = "harmony", dims = 1:50)
SO.harmony <- FindNeighbors(SO.harmony, dims = 1:50, reduction = "harmony")

SO.harmony <- FindClusters(SO.harmony, resolution = 1, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.9, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.8, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.7, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.6, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.5, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.4, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.3, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.2, reduction = "harmony")
SO.harmony <- FindClusters(SO.harmony, resolution = 0.1, reduction = "harmony")

setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_Stuff/CD8")
save(SO.harmony, file = "CD8_harmony_OnlyTCRscells_18062024.RData")

###SELECTING RESOLUTION TO FURTHER ANALYZE#### FIGURES 1 f,g,h Author: Emilia Schneider Revueltas

Idents(SO.harmony) <- "RNA_snn_res.0.7"

UMAP <- DimPlot(SO.harmony, pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP, id = "ident", size = 4, box = T, repel = F) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP2 <- DimPlot(SO.harmony, split.by = "donor", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP2, id = "ident", box = T, size = 3) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP3 <- DimPlot(SO.harmony, split.by = "celltype", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP3, id = "ident", box = T, size = 3) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

U2_U3 <- (UMAP2 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))) / (UMAP3 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)))
U2_U3 
