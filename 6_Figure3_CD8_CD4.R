library(SeuratObject)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)



#For Fig 3, change pathfile and CD4 and CD8 annotations

####VOLCANOS CLUSTER BY CLUSTER####
#Bind resolution and celltypes CD8
SO@meta.data$res7 = SO@meta.data$RNA_snn_res.0.7
table(SO@meta.data$res7, SO@meta.data$celltype)

SO@meta.data$fn_celltype <- as.factor(paste(SO@meta.data$res7, SO@meta.data$celltype, sep = "_"))
colnames(SO@meta.data)
class(SO@meta.data$fn_celltype)

Idents(SO) <- "fn_celltype"

#Input the clusters to be compared
ID1 = "1_Blood_CD8p_CD69n"
ID2 = "1_BM_CD8p_CD69n"

ID1_ID2 <- FindMarkers(SO, ident.1 = ID1, 
                       ident.2 = ID2, only.pos = F, min.pct = 0.3)
ID1_ID2$gene <- rownames(ID1_ID2)
volcano<-ID1_ID2

file_path = paste("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8/Paper/Figure_3/Volcano_", ID1, "_vs_", ID2, ".csv", sep = "")
print(file_path)
write.csv(ID1_ID2, 
          file = file_path)

myVolcanoPlot(volcano, paste(ID1, "vs.", ID2))

#Function for Volcano Plot
myVolcanoPlot <- function(volcano, label) {
  
  volcano$gene <- rownames(volcano)
  volcano$diffexp <- "NO"
  volcano$diffexp[volcano$avg_log2FC > +0.35 & volcano$p_val < 0.05] <- "UP"
  volcano$diffexp[volcano$avg_log2FC < -0.35 & volcano$p_val < 0.05] <- "DOWN"
  volcano$plotlabel <- NA
  volcano$plotlabel[volcano$diffexp != "NO"] <- volcano$gene[volcano$diffexp != "NO"]
  volcano$plotlabel[grep("^RP", x=rownames(volcano))] <- NA
  volcano$plotlabel[grep("^MT", x=rownames(volcano))] <- NA
  
  ggplot(data=volcano, aes(x=avg_log2FC, y= -log10(p_val), 
                           col = diffexp, label = plotlabel)) + 
    geom_point() + 
    theme_minimal() + 
    geom_vline(xintercept=c(-0.35, 0.35), col="grey", linetype = "dashed") +
    #geom_hline(yintercept=-log10(0.05), col="red") + 
    geom_text_repel(max.overlaps = 15, size = 7) +
    xlab(label = label) +
    scale_y_continuous(limits = c(0, 51)) +
    scale_x_continuous(limits = c(-2.5, 1)) +
    theme_classic(base_size = 25) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +  
    NoLegend()
}


##Bubble plot per tissue##

Idents(SO) <- "celltype"

#MIN PCT 05 and LFC 0.25
CD8h <- FindAllMarkers(SO, only.pos = TRUE, 
                       min.pct = 0.5, 
                       logfc.threshold = 0.25)
write.csv(CD8h, 
          file ="/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8/Paper/Figure_3/CD8_res7_celltype_FAM_OTCR_opT_5_25.csv")


CD8h_FAM_top30 <- CD8h %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
CD8h_FAM_top30 <- as.factor(unique(CD8h_FAM_top30$gene))
features = CD8h_FAM_top30

DotPlot(SO, features = features, dot.min = 0) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  scale_color_gradientn(colours = c("#7066bc","white","#56ae6c"))




