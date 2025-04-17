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
library(eulerr)
library(reshape2)

#install.packages("harmony")


setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4")


file_name = "CD4_harmony_OnlyTCRscells_18062024.RData"


load(file_name)

SO <- SO.harmony

setwd("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4/Paper/Res6")

####Pawel Freddy Colors ####

Disc_colors = c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
                "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",
                "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")

Disc_colors_cust = c("#999999","#999999","#1965B0","#999999","#882E72","#999999","#999999","#FDB462","#999999","#999999",
                     "#999999","#999999")

Disc_colors_cust2 = c("#999999",  "#FB8072", "#999999", "#999999", "#999999","#B17BA6","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999","#999999",
                      "#999999","#999999")

#### Choose IDENTS####
#Run only desired IDENT
Idents(SO) <- "RNA_snn_res.0.6"
DefaultAssay(SO) <- "RNA"

# # # # # # # # # # # # # # # # # # # #

####Basic UMAP#### Figures 1 a,b,c

# # # # # # # # # # # # # # # # # # # #

UMAP <- DimPlot(SO, pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP, id = "ident", size = 4, box = T, repel = F) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP <- DimPlot(SO, pt.size = 0.3, cols = Disc_colors_cust2)
LabelClusters(UMAP, id = "ident", size = 4, box = T, repel = F) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP2 <- DimPlot(SO, split.by = "donor", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP2, id = "ident", size = 3, box = T, repel = F) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

UMAP3 <- DimPlot(SO, split.by = "celltype", pt.size = 0.3, cols = Disc_colors)
LabelClusters(UMAP3, id = "ident", size = 3, box = T, repel = F) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))

U2_U3 <- (UMAP2 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2))) / (UMAP3 + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)))
U2_U3 

# # # # # # # # # # # # # # # # # # # #

####Figure 1d,e###

# # # # # # # # # # # # # # # # # # # #

####Cluster contribution per DONOR#### Fig 1d

#BAR PLOT CELL NUMBERS PER CLUSTER OF EACH DONOR
Celln <- table(SO@meta.data$RNA_snn_res.0.6, SO@meta.data$donor)
Celln_df <- as.data.frame.matrix(Celln)
cluster_values <- 0:11

Celln_df$cluster_values <- rownames(Celln_df)
Celln_df$cluster_values <- factor(Celln_df$cluster_values, levels = as.character(0:11))
print(Celln_df)

# Reshape the data frame to long format for ggplot2
Celln_long <- Celln_df %>%
  gather(key = "Donor", value = "Count", -cluster_values)

# Create the bar plot with absolute cell counts
ggplot(Celln_long, aes(x = cluster_values, y = Count, fill = Donor)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = Disc_colors[5:7]) +
  labs(title = 'Counts by Cluster and Donor', x = 'Cluster', y = 'Count') +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )

####Percentages

Celln_percent_df <- Celln_df %>%
  mutate(across(D1:D3, ~ . / sum(.))) %>%
  mutate(across(D1:D3, ~ . * 100))

# Convert cluster_values to a factor with levels in the correct order
Celln_percent_df$cluster_values <- factor(Celln_percent_df$cluster_values, levels = as.character(0:11))

# View the resulting data frame
print(Celln_percent_df)

# Reshape the data frame to long format for ggplot2
Celln_long <- Celln_percent_df %>%
  gather(key = "Donor", value = "Percentage", -cluster_values)


mean_data <- Celln_long %>%
  group_by(cluster_values) %>%
  summarise(mean_Percentage = mean(Percentage, na.rm = TRUE))

###DOTPLOT with a bar for representing the mean
ggplot() +
  geom_point(data = Celln_long, aes(x = cluster_values, y = Percentage, colour = Donor),
             stat = 'identity', position = position_dodge(width = 0.2), size = 3) +  # Points
  geom_bar(data = mean_data, aes(x = cluster_values, y = mean_Percentage),
           stat = "identity", fill = "grey", alpha = 0.5, width = 0.6, position = "dodge") +  # Mean bars
  scale_color_manual(values = Disc_colors[5:7]) +
  labs(title = 'Cluster contribution per Donor', x = 'Cluster', y = 'Frequencies among donor-s total cells') +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )

###DOTPLOT WITH BARS FOR THE MEAN
ggplot() +
  geom_point(data = Celln_long, aes(x = cluster_values, y = Percentage, colour = Donor),
             stat = 'identity', position = position_dodge(width = 0.2), size = 3) +  # Points
  geom_segment(data = mean_data, aes(x = cluster_values, xend = cluster_values,
                                     y = mean_Percentage - 0.1, yend = mean_Percentage + 0.1),
               colour = "black", linewidth = 7.7) +  # Mean lines
  scale_color_manual(values = Disc_colors[5:7]) +
  labs(title = 'Cluster contribution per Donor', x = 'Cluster', y = 'Frequencies among donor-s total cells') +
  coord_cartesian(ylim = c(min(Celln_long$Percentage, na.rm = TRUE) - 0.1, 
                           max(Celln_long$Percentage, na.rm = TRUE) + 0.1)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )

####Cluster contribution per CELLTYPE### Fig 1 e

#BAR PLOT CELL NUMBERS PER CLUSTER OF EACH celltype
Celln_type <- table(SO@meta.data$RNA_snn_res.0.6, SO@meta.data$celltype)
Celln_type_df <- as.data.frame.matrix(Celln_type)

Celln_type_df$cluster_values <- rownames(Celln_type_df)
Celln_type_df$cluster_values <- factor(Celln_type_df$cluster_values, levels = as.character(0:11))
print(Celln_type_df)

# Reshape the data frame to long format for ggplot2
Celln_type_long <- Celln_type_df %>%
  gather(key = "Origin", value = "Count", -cluster_values)

# Create the bar plot with absolute cell counts
ggplot(Celln_type_long, aes(x = cluster_values, y = Count, fill = Origin)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = Disc_colors[5:7]) +
  labs(title = 'Counts by Cluster and Donor', x = 'Cluster', y = 'Count') +
  theme_minimal()+
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )

####Percentages

Celln_type_percent_df <- Celln_type_df %>%
  mutate(across(Blood_CD4p_CD69n:BM_CD4p_CD69p, ~ . / sum(.))) %>%
  mutate(across(Blood_CD4p_CD69n:BM_CD4p_CD69p, ~ . * 100))

# Convert cluster_values to a factor with levels in the correct order
Celln_type_percent_df$cluster_values <- factor(Celln_type_percent_df$cluster_values, levels = as.character(0:11))

# View the resulting data frame
print(Celln_type_percent_df)

# Reshape the data frame to long format for ggplot2
Celln_type_long <- Celln_type_percent_df %>%
  gather(key = "Donor", value = "Percentage", -cluster_values)


mean_data <- Celln_type_long %>%
  group_by(cluster_values) %>%
  summarise(mean_Percentage = mean(Percentage, na.rm = TRUE))

###DOTPLOT with a bar for representing the mean
ggplot() +
  geom_point(data = Celln_type_long, aes(x = cluster_values, y = Percentage, colour = Donor),
             stat = 'identity', position = position_dodge(width = 0.2), size = 3) +  # Points
  #geom_bar(data = mean_data, aes(x = cluster_values, y = mean_Percentage),
  #stat = "identity", fill = "grey", alpha = 0.5, width = 0.6, position = "dodge") +  # Mean bars
  scale_color_manual(values = Disc_colors[9:11]) +
  labs(title = 'Cluster contribution per Origin', x = 'Cluster', y = 'Frequencies among origin total cells') +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )


# # # # # # # # # # # # # # # # # # # #

####FIND ALL MARKERS#### Figure 1f

# # # # # # # # # # # # # # # # # # # #


DefaultAssay(SO.harmony) <- "RNA"
Idents(SO.harmony) <- "RNA_snn_res.0.6"
#MIN PCT 05 and LFC 0.25
CD4h_FAM <- FindAllMarkers(SO.harmony, only.pos = TRUE, 
                           min.pct = 0.5, 
                           logfc.threshold = 0.25)
write.csv(CD4h_FAM, 
          file ="/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4/Paper/Res6/CD4_res6_FAM_OTCR_opT_5_25.csv")

CD4h_FAM_top10 <- CD4h_FAM %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
CD4h_FAM_top10 <- as.factor(unique(CD4h_FAM_top10$gene))
CD4h_FAM_top7 <- CD4h_FAM %>% group_by(cluster) %>% top_n(n = 7, wt = avg_log2FC)
CD4h_FAM_top7 <- as.factor(unique(CD4h_FAM_top7$gene))
CD4h_FAM_top5 <- CD4h_FAM %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
CD4h_FAM_top5 <- as.factor(unique(CD4h_FAM_top5$gene))


#MIN PCT 01 and LFC 0.25
CD4h_FAM_1_3 <- FindAllMarkers(SO.harmony, only.pos = TRUE, 
                               min.pct = 0.1, 
                               logfc.threshold = 0.3)

write.csv(CD4h_FAM_1_3, 
          file ="/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4/Paper/Res6/CD4_res6_FAM_OTCR_opT_1_3.csv")

CD4h_FAM_1_3_top10 <- CD4h_FAM_1_3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
CD4h_FAM_1_3_top10 <- as.factor(unique(CD4h_FAM_1_3_top10$gene))


#PLOT: introduce the vector of genes to be plotted
features = CD4h_FAM_1_3_top10
features = CD4h_FAM_top10
features = FEATURES_CD4
features = CD4h_FAM_top7

DotPlot(SO.harmony, features = features, dot.min = 0) +
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  scale_color_gradientn(colours = c("#FFC000","white","#F50087"))
