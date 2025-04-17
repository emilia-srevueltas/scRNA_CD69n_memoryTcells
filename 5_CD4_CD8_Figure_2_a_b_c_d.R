#Resolution has to be changed accordingly


####SHANNON INDEX FUNCTION####
Disc_colors_1 = c("#DC050C", "#FB8072", "#33A02C","#B2DF8A", "#1965B0", "#7BAFDE", "#882E72","#B17BA6","#FF7F00","#FDB462"
                  ,"#E7298A","#E78AC3","#55A1B1","#8DD3C7","#A6761D","#E6AB02","#7570B3","#BEAED4","#666666","#999999",
                  "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000","#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")



GetPOP <- function(donor, cluster) { 
  donor_filter = SO@meta.data$donor == donor
  cluster_filter = SO@meta.data$RNA_snn_res.0.6[donor_filter] == cluster
  
  Clones <- as.factor(SO@meta.data$aa_clones[donor_filter][cluster_filter])
  freq_table <- table(Clones)
  freq_table = freq_table[freq_table > 0]
  sorted_freq <- sort(freq_table, decreasing = TRUE)
  sorted_freq_df <- as.data.frame(sorted_freq)
  pop = sorted_freq_df$Freq
  return(pop)
}

ShannonIndex <- function(pop) { 
  n = pop[pop>0]
  N = sum(pop)
  p = n/N
  H = -sum(p*log(p))
  return(H)
}




# Define cluster values
cluster_values <- 0:11

# Define donor values
donors <- c("D1", "D2", "D3")

# Initialize an empty data frame to store the results
Shannon_Indexes <- data.frame(matrix(nrow = length(cluster_values), ncol = length(donors)))
colnames(Shannon_Indexes) <- donors
rownames(Shannon_Indexes) <- cluster_values

# Loop through each donor
for (donor in donors) {
  # Initialize an empty list to store Shannon Index values for the current donor
  ShannonList <- c()
  
  # Loop through each cluster value for the current donor
  for (cluster in cluster_values) {
    # Compute Shannon Index for the current donor and cluster
    Shannon_value <- ShannonIndex(GetPOP(donor, cluster))
    
    # Append the computed Shannon Index value to ShannonList
    ShannonList <- c(ShannonList, Shannon_value)
  }
  
  # Add ShannonList as a column to Shannon_Indexes with the donor name as the column name
  Shannon_Indexes[, donor] <- ShannonList
}

write.csv(Shannon_Indexes, 
          file = "/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4/Paper/Res6/Shannon_Indexes_res6.csv")

###ShannonBarPlot
Shannon_Indexes_df <- data.frame(
  cluster_values = cluster_values,
  D1 = Shannon_Indexes$D1,
  D2 = Shannon_Indexes$D2,
  D3 = Shannon_Indexes$D3)

# Plot the data
Shannon_Indexes_long <- pivot_longer(Shannon_Indexes_df, cols = starts_with("D"), names_to = "donor", values_to = "ShannonList")

# Create the stacked bar plot
ggplot(Shannon_Indexes_long, aes(x = as.factor(cluster_values), y = ShannonList, fill = donor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = Disc_colors[5:7]) +
  labs(title = "Shannon Entropy of each Cluster", x = "Cluster", y = "Shannon Entropy") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )


####Simpson Index FUNCTION####
GetPOP <- function(donor, cluster) { 
  donor_filter = SO@meta.data$donor == donor
  cluster_filter = SO@meta.data$RNA_snn_res.0.6[donor_filter] == cluster
  
  Clones <- as.factor(SO@meta.data$aa_clones[donor_filter][cluster_filter])
  freq_table <- table(Clones)
  freq_table = freq_table[freq_table > 0]
  sorted_freq <- sort(freq_table, decreasing = TRUE)
  sorted_freq_df <- as.data.frame(sorted_freq)
  pop = sorted_freq_df$Freq
  return(pop)
}


SimpsonIndex <- function(pop) {
  
  p <- pop / sum(pop)
  # Compute the Simpson Index
  D <- 1 - sum(p^2)
  return(D)
}


#Define
cluster_values <- 0:11
donors <- c("D1", "D2", "D3")

# Initialize an empty data frame with the appropriate number of rows
Simpson_Indexes <- data.frame(matrix(nrow = length(cluster_values), ncol = length(donors)))
colnames(Simpson_Indexes) <- donors
rownames(Simpson_Indexes) <- cluster_values

# Loop through each donor
for (donor in donors) {
  # Initialize an empty list to store Simpson Index values for the current donor
  SimpsonList <- c()
  
  # Loop through each cluster value for the current donor
  for (cluster in cluster_values) {
    # Compute Simpson Index for the current donor and cluster
    Simpson_value <- SimpsonIndex(GetPOP(donor, cluster))
    
    # Append the computed Simpson Index value to SimpsonList
    SimpsonList <- c(SimpsonList, Simpson_value)
  }
  
  # Add SimpsonList as a column to Simpson_Indexes with the donor name as the column name
  Simpson_Indexes[, donor] <- SimpsonList
}

write.csv(Simpson_Indexes, 
          file = "/Users/emiliasere/Desktop/Emilia/scRNA/R_stuff/CD4/Paper/Res6/Simpson_Indexes_res6.csv")

###pLOT

Simpson_Indexes_df <- data.frame(
  cluster_values = cluster_values,
  D1 = Simpson_Indexes$D1,
  D2 = Simpson_Indexes$D2,
  D3 = Simpson_Indexes$D3)

# Plot the data
Simpson_Indexes_long <- pivot_longer(Simpson_Indexes_df, cols = starts_with("D"), names_to = "donor", values_to = "SimpsonList")

# Create the stacked bar plot
ggplot(Simpson_Indexes_long, aes(x = as.factor(cluster_values), y = SimpsonList, fill = donor)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = Disc_colors[5:7]) +
  labs(title = "Simpson of each Cluster", x = "Cluster", y = "Simpson Index") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, hjust = 1),
    axis.text.x = element_text(size = 18, hjust = 1),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_blank()  # Remove legend title if not needed
  )

