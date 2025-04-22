library(ggplot2)
library(patchwork)
library(SeuratObject)
library(Seurat)

##### FUNCTION DEFINITIONS #####

GetOverlapBetweenTwo <- function(labels, groups, group_1, group_2)
{
  temp_labels <- labels[groups == group_1 | groups == group_2]
  temp_groups <- groups[groups == group_1 | groups == group_2]
  
  group_1 <- which(temp_groups == group_1)
  group_2 <- which(temp_groups == group_2)
  
  real_overlap <- GetOverlapOfTwo(temp_labels[group_1], temp_labels[group_2])
  
  shuffled_overlaps <- c()
  for (i in c(1:1000))
  {
    temp_labels <- sample(temp_labels)
    shuffled_overlap <- GetOverlapOfTwo(temp_labels[group_1], temp_labels[group_2])
    
    shuffled_overlaps <- c(shuffled_overlaps, shuffled_overlap)
  }
  
  return(list(real_overlap=real_overlap, shuffled_overlaps=shuffled_overlaps))
}

GetOverlapBetweenThree <- function(labels, groups, group_1, group_2, group_3)
{
  temp_labels <- labels[groups == group_1 | groups == group_2 | groups == group_3]
  temp_groups <- groups[groups == group_1 | groups == group_2 | groups == group_3]
  
  group_1 <- which(temp_groups == group_1)
  group_2 <- which(temp_groups == group_2)
  group_3 <- which(temp_groups == group_3)
  
  real_overlap <- GetOverlapOfThree(temp_labels[group_1], temp_labels[group_2], temp_labels[group_3])
  
  shuffled_overlaps <- c()
  for (i in c(1:1000))
  {
    temp_labels <- sample(temp_labels)
    shuffled_overlap <- GetOverlapOfThree(temp_labels[group_1], temp_labels[group_2], temp_labels[group_3])
    shuffled_overlaps <- c(shuffled_overlaps, shuffled_overlap)
  }
  
  return(list(real_overlap=real_overlap, shuffled_overlaps=shuffled_overlaps))
}

GetOverlapOfTwo <- function(elements_1, elements_2)
{
  return(length(intersect(unique(elements_1), unique(elements_2))))
}

GetOverlapOfThree <- function(elements_1, elements_2, elements_3)
{
  return(length(intersect(intersect(unique(elements_1), unique(elements_2)), unique(elements_3))))
}

GetPVals <- function(overlaps)
{
  real <- overlaps$real_overlap
  all <- overlaps$shuffled_overlaps
  N <- length(all)
  
  p_lt <- length(which(all <= real)) / N
  p_gt <- length(which(all >= real)) / N
  
  return(list(p_lt = p_lt, p_gt = p_gt))
}

##### OVERLAP GENERATION #####

clonotypes <- SO$aa_clones
celltypes <- SO$celltype
donors <- SO$donor
clusters <- SO$RNA_snn_res.0.7

#donor_names <- c("D1") 
#cluster_names <- c("4","5")

donor_names <- c("D1", "D2", "D3")
cluster_names <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")


for(donor in donor_names)
{
  Overlap_Matrix <- data.frame(matrix(ncol = length(cluster_names), nrow = 12))
  rownames(Overlap_Matrix) <- c("overlap_Bl_BMn", "p_lt_Bl_BMn", "p_gt_Bl_BMn", 
                                "overlap_Bl_BMp", "p_lt_Bl_BMp", "p_gt_Bl_BMp",
                                "overlap_BMn_BMp", "p_lt_BMn_BMp", "p_gt_BMn_BMp",
                                "overlap_Bl_BMn_BMp", "p_lt_Bl_BMn_BMp", "p_gt_Bl_BMn_BMp")
  colnames(Overlap_Matrix) <- cluster_names
  
  
  Cluster_ovs <- c()
  for(cluster in cluster_names)
  {
    
    filter <- which(donors == donor & clusters == cluster)
    
    clonotypes_temp <- clonotypes[filter]
    celltypes_temp <- celltypes[filter]
    
    overlap_Bl_BMn <- GetOverlapBetweenTwo(clonotypes_temp, celltypes_temp, "Blood_CD8p_CD69n", "BM_CD8p_CD69n")
    overlap_Bl_BMp <- GetOverlapBetweenTwo(clonotypes_temp, celltypes_temp, "Blood_CD8p_CD69n", "BM_CD8p_CD69p")
    overlap_BMn_BMp <- GetOverlapBetweenTwo(clonotypes_temp, celltypes_temp, "BM_CD8p_CD69n", "BM_CD8p_CD69p")
    overlap_Bl_BMn_BMp <- GetOverlapBetweenThree(clonotypes_temp, celltypes_temp, "Blood_CD8p_CD69n", "BM_CD8p_CD69n", "BM_CD8p_CD69p")
    
    p_vals_Bl_BMn <- GetPVals(overlap_Bl_BMn)
    p_vals_Bl_BMp <- GetPVals(overlap_Bl_BMp)
    p_vals_BMn_BMp <- GetPVals(overlap_BMn_BMp)
    p_vals_Bl_BMn_BMp <- GetPVals(overlap_Bl_BMn_BMp)
    
    Cluster_ov <- c(overlap_Bl_BMn[["real_overlap"]], p_vals_Bl_BMn[["p_lt"]], p_vals_Bl_BMn[["p_gt"]],
                    overlap_Bl_BMp[["real_overlap"]], p_vals_Bl_BMp[["p_lt"]], p_vals_Bl_BMp[["p_gt"]],
                    overlap_BMn_BMp[["real_overlap"]], p_vals_BMn_BMp[["p_lt"]], p_vals_BMn_BMp[["p_gt"]],
                    overlap_Bl_BMn_BMp[["real_overlap"]], p_vals_Bl_BMn_BMp[["p_lt"]], p_vals_Bl_BMn_BMp[["p_gt"]])
    
    Cluster_ovs <- c(Cluster_ovs, Cluster_ov)
    
    p <- ggplot(data.frame(overlap_Bl_BMn[["shuffled_overlaps"]]), aes(x = overlap_Bl_BMn[["shuffled_overlaps"]])) +
      scale_y_continuous(expand = c(0, 0), name = "counts") +
      geom_histogram(binwidth = 1, fill = "#E78AC3", color = "#e9ecef", alpha = 0.9) +
      ggtitle(paste(donor, cluster,"1000 randomizations \n BLCD8+CD69- BMCD8+CD69-")) +
      theme(plot.title = element_text(size = 15)) +
      geom_vline(xintercept = overlap_Bl_BMn[["real_overlap"]], col = "black") +
      annotate("text", x = overlap_Bl_BMn[["real_overlap"]], y = 0, 
               label = overlap_Bl_BMn[["real_overlap"]], vjust = 0, col = "red")+
      annotate("text", x = Inf, y = Inf, label = paste("p = ", p_vals_Bl_BMn[["p_lt"]]), 
               hjust = 1, vjust = 1, col = "black", size = 4, fontface = "bold")
    
    p1 <- ggplot(data.frame(overlap_Bl_BMp[["shuffled_overlaps"]]), aes(x = overlap_Bl_BMp[["shuffled_overlaps"]])) +
      scale_y_continuous(expand = c(0, 0), name = "counts") +
      geom_histogram(binwidth = 1, fill = "#E78AC3", color = "#e9ecef", alpha = 0.9) +
      ggtitle(paste(donor, cluster,"1000 randomizations \n BLCD8+CD69- BMCD8+CD69+")) +
      theme(plot.title = element_text(size = 15)) +
      geom_vline(xintercept = overlap_Bl_BMp[["real_overlap"]], col = "black") +
      annotate("text", x = overlap_Bl_BMp[["real_overlap"]], y = 0, 
               label = overlap_Bl_BMp[["real_overlap"]], vjust = 0, col = "red")+
      annotate("text", x = Inf, y = Inf, label = paste("p = ", p_vals_Bl_BMp[["p_lt"]]), 
               hjust = 1, vjust = 1, col = "black", size = 4, fontface = "bold")
    
    p2 <- ggplot(data.frame(overlap_BMn_BMp[["shuffled_overlaps"]]), aes(x = overlap_BMn_BMp[["shuffled_overlaps"]])) +
      scale_y_continuous(expand = c(0, 0), name = "counts") +
      geom_histogram(binwidth = 1, fill = "#E78AC3", color = "#e9ecef", alpha = 0.9) +
      ggtitle(paste(donor, cluster,"1000 randomizations \n BMCD8+CD69- BMCD8+CD69+")) +
      theme(plot.title = element_text(size = 15)) +
      geom_vline(xintercept = overlap_BMn_BMp[["real_overlap"]], col = "black") +
      annotate("text", x = overlap_BMn_BMp[["real_overlap"]], y = 0, 
               label = overlap_BMn_BMp[["real_overlap"]], vjust = 0, col = "red")+
      annotate("text", x = Inf, y = Inf, label = paste("p = ", p_vals_BMn_BMp[["p_lt"]]), 
               hjust = 1, vjust = 1, col = "black", size = 4, fontface = "bold")
    
    p3 <- ggplot(data.frame(overlap_Bl_BMn_BMp[["shuffled_overlaps"]]), aes(x = overlap_Bl_BMn_BMp[["shuffled_overlaps"]])) +
      scale_y_continuous(expand = c(0, 0), name = "counts") +
      geom_histogram(binwidth = 1, fill = "#E78AC3", color = "white", alpha = 0.9) +
      ggtitle(paste(donor, cluster,"1000 randomizations \n BLCD8+CD69- BMCD8+CD69- BMCD8+CD69+")) +
      theme(plot.title = element_text(size = 15)) +
      geom_vline(xintercept = overlap_Bl_BMn_BMp[["real_overlap"]], col = "black") +
      annotate("text", x = overlap_Bl_BMn_BMp[["real_overlap"]], y = 0, 
               label = overlap_Bl_BMn_BMp[["real_overlap"]], vjust = 0, col = "red")+
      annotate("text", x = Inf, y = Inf, label = paste("p = ", p_vals_Bl_BMn_BMp[["p_lt"]]),
               #"\nOv:",overlap_Bl_BMn_BMp[["real_overlap"]]), 
               hjust = 1, vjust = 1, col = "black", size = 4, fontface = "bold")
    
    PLOT <- (p | p1) / (p2 | p3)
    
    pdf_filename <- paste0(donor, "_cluster", cluster, "_overlap.pdf")
    ggsave(pdf_filename, PLOT, device = "pdf", width = 11.69, height = 8.27, dpi = 300)
  }
  
  Overlap_Matrix[,cluster_names] <- c(Cluster_ovs)
  
  file_path = paste("/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD8/Paper/Res7/Overlapps_12092024/Ov_Matrix_", donor, ".csv", sep = "")
  print(file_path)
  write.csv(Overlap_Matrix, file = file_path)
  
}

