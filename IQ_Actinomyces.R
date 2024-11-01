
# Packages instalation

if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("RColorBrewer","tidyverse", "readxl", "rvest","dplyr","tidyr","igraph","visNetwork", "readr","ggplot2", "ggraph", "graphTweets", "writexl")
#install.packages("rJava")
#install.packages("rcdk")
#install.packages("proxy")
Sys.getenv("JAVA_HOME")
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk-17')
library(rJava)
library(rcdk)
library(httr)

# List of directories to create
Dirs <- c("Data/CSV_IQ", "Data/IQ_Plot","Data/CSV_IQNPA", "Data/IQNPA_Plot")

# Loop through the list and create directories if they don't exist
for (dir in Dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Load network from a GraphML file. Change the file name as necessary
graphml <- read_graph("Actino_network.graphml", format = "graphml")
Metadata <- read_delim("Metadata_Actinomycetes.csv", 
                       delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Get node and edge data from the graph
if (is_igraph(graphml)) {
  df_edges <- as_data_frame(graphml, what = "edges")
  df_nodes <- as_data_frame(graphml, what = "vertices")
}

## Calculate total nodes per sample
Total_nodes_per_strain <- df_nodes %>%
  select(id, ATTRIBUTE_Strain, MQScore) %>% 
  separate_rows(ATTRIBUTE_Strain, sep = "\\,") %>%  
  group_by(ATTRIBUTE_Strain)

# Count total nodes per strain
Total_nodes <- table(Total_nodes_per_strain$ATTRIBUTE_Strain)
Total_nodes <- as.data.frame(Total_nodes)
colnames(Total_nodes) <- c("ATTRIBUTE_strain", "Count")

# Extract index inputs for MQS
score_MQS <- Total_nodes_per_strain %>%
  filter(!is.na(MQScore) & MQScore != "")

# Count nodes with annotated strain information
nodes_annotated_strain <- table(score_MQS$ATTRIBUTE_Strain)
nodes_annotated_strain <- as.data.frame(nodes_annotated_strain)
colnames(nodes_annotated_strain) <- c("ATTRIBUTE_strain", "Count")

#Jump to line 910 for IQ based on NPAtlas anottations


## Define strata for MQScore
MQS_breaks <- seq(0.00, 1, by = 0.05)
strata_MQS <- setNames(
  lapply(seq_along(MQS_breaks)[-length(MQS_breaks)], function(i) c(MQS_breaks[i], MQS_breaks[i + 1])),
  paste0(sprintf("%.2f", MQS_breaks[-length(MQS_breaks)]), "-", sprintf("%.2f", MQS_breaks[-1])))

## Function to stratify and count values, ignoring NAs
stratify_and_count_unique <- function(df, column, strata) {
  unique_df <- df %>% 
    group_by(id) %>% 
    slice(1) %>%  # Count only the first occurrence of each unique identifier (id) in the specified column
    ungroup()
  
  # Count unique values for each stratum
  counts <- sapply(strata, function(range) {
    sum(unique_df[[column]] >= range[1] & unique_df[[column]] < range[2], na.rm = TRUE)
  })
  
  return(counts)
}

## Create lists to store individual dataframes and strata counts
individual_dfs <- list()
strata_counts <- list()

## For each sample, create an individual dataframe and calculate strata counts
for (sample_name in unique(score_MQS$ATTRIBUTE_Strain)) {
  individual_df <- score_MQS %>%
    filter(ATTRIBUTE_Strain == sample_name) %>%
    select(id, MQScore)
  
  # Count unique values in each stratum for MQScore
  MQScore_counts <- stratify_and_count_unique(individual_df, "MQScore", strata_MQS)
  
  # Ensure that the strata names are consistent
  counts_df <- data.frame(
    Stratum = names(strata_MQS),
    MQScore_Count = MQScore_counts
  )
  
  # Store the count dataframe in the list
  strata_counts[[sample_name]] <- counts_df
  
  # Assign the individual dataframe to the list
  individual_dfs[[sample_name]] <- individual_df
}

## Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder
for(sample_name in names(individual_dfs)) {
  write.csv(individual_dfs[[sample_name]], file = paste0("Data/CSV_IQ/", sample_name, "_data_MQS.csv"), row.names = FALSE)
  
  
  # Assign each dataframe to a variable in the R environment
  assign(paste0(sample_name, "_data"), individual_df)
  assign(paste0(sample_name, "_strata_counts"), counts_df)
}

## Create a data frame to store strain names
strains <- data.frame(Total_nodes$ATTRIBUTE_strain)
colnames(strains) <- c("strain")
strains <- strains$strain

## Create an empty dataframe for the summary
shannon_summary_MQS <- data.frame(Strain = character(), Shannon_Index = numeric(), stringsAsFactors = FALSE)

## Iterate over each strain
for (strain in strains) {
  # Name of the corresponding dataframe
  df_name <- paste0(strain, "_strata_counts")
  
  # Check if the dataframe exists
  if (exists(df_name)) {
    df_strata_counts <- get(df_name)
    
    # Obtain the total count for this strain from Total_nodes
    total <- Total_nodes$Count[Total_nodes$ATTRIBUTE_strain == strain]
    
    # Check that the total is greater than zero
    if (total > 0) {
      # Calculate the proportion (pi) for each stratum and add it as a new column
      df_strata_counts$pi <- df_strata_counts$MQScore_Count / total
      
      # Add a column with the total count
      df_strata_counts$Total <- total
      
      # Calculate pi * ln(pi), handling the case where pi = 0
      df_strata_counts$pi_ln_pi <- ifelse(df_strata_counts$pi > 0, df_strata_counts$pi * log(df_strata_counts$pi), 0)
    } else {
      # If the total is zero, assign NA to the new columns
      df_strata_counts$pi <- NA
      df_strata_counts$Total <- NA
      df_strata_counts$pi_ln_pi <- NA
    }
    
    # Calculate the Shannon index as -sum(pi * ln(pi))
    shannon_index <- -sum(df_strata_counts$pi_ln_pi, na.rm = TRUE)
    
    # Add the Shannon index as a new column (same for all strata)
    df_strata_counts$shannon_index <- shannon_index
    
    # Append the Shannon index to the summary
    shannon_summary_MQS <- rbind(shannon_summary_MQS, data.frame(Strain = strain, Shannon_Index = shannon_index))
    
    # Select only the relevant columns before saving
    df_strata_counts <- df_strata_counts[, c("Stratum", "MQScore_Count", "pi", "Total", "pi_ln_pi", "shannon_index")]
    
    # Save the updated dataframe to a CSV file
    write.csv(df_strata_counts, file = paste0("Data/CSV_IQ/", df_name, "_SI_MQS.csv"), row.names = FALSE)
  } else {
    warning(paste("The dataframe", df_name, "does not exist."))
  }
}


## Save Shannon Index summary to a CSV file
write.csv(shannon_summary_MQS, file = "Data/CSV_IQ/1_MQS_shannon_index_summary.csv", row.names = FALSE)

# Extract index input: COS_score

# Get cosine scores where node2 matches the id in score_MQS
score_COS_node1 <- df_edges %>%
  select(node1, cosine_score, node2) %>%
  semi_join(score_MQS, by = c("node2" = "id")) # Retrieve cosine values from node to node

# Get complementary cosine scores where node1 matches the id in score_MQS
score_COS_node2 <- df_edges %>%
  select(node2, cosine_score, node1) %>%
  semi_join(score_MQS, by = c("node1" = "id")) # Retrieve complementary cosine values

## Combine interactions into a single dataframe
combined_interactions <- bind_rows(
  score_COS_node1 %>% mutate(direction = "Score_COS_node1"),
  score_COS_node2 %>% mutate(direction = "Score_COS_node2")
) # Add a new column called direction to indicate the source node of the interaction.

# Create a combined node column based on direction
combined_interactions <- combined_interactions %>% mutate(
  combined_nodes = ifelse(direction == "Score_COS_node1",
                          paste(node1, node2, sep = "/"),
                          paste(node2, node1, sep = "/"))
) # If “Score_COS_node1”, concatenate node1 and node2; else, concatenate node2 and node1.

# Separate combined nodes into nodeA and nodeB
combined_interactions <- combined_interactions %>%
  separate(combined_nodes, into = c("nodeA", "nodeB"), sep = "/") # Assign nodeA as the node of interest.

## Remove rows where nodeA and nodeB are the same (singlets)
combined_interactions <- combined_interactions[combined_interactions$nodeA != combined_interactions$nodeB, ]

# Select relevant columns for cosine scores
score_COS <- combined_interactions %>%
  select(nodeA, cosine_score) # Filter cosine values for the node of interest.

# Merge cosine scores with node metadata
score_COS <- merge(score_COS, df_nodes, by.x = "nodeA", by.y = "id", all = TRUE)  
# Note: MQS values are repeated when there are multiple COS scores, as intended.

# Clean up score_COS dataframe
score_COS <- score_COS %>% select(nodeA, cosine_score, ATTRIBUTE_Strain) %>% 
  filter(!is.na(cosine_score) & cosine_score != "") %>% 
  separate_rows(ATTRIBUTE_Strain, sep = "\\,") %>%  
  group_by(ATTRIBUTE_Strain)

# Rename columns for clarity
colnames(score_COS) <- c("id", "Cosine_score", "ATTRIBUTE_Strain")

# Merge with metadata to include additional information
score_COS <- merge(score_COS, Metadata, 
                   by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain", all.x = TRUE)

# Final selection of relevant columns
score_COS <- score_COS %>% select(id, Cosine_score, ATTRIBUTE_Strain)

## Define strata for cosine score
cos_breaks <- seq(0.70, 1.00, by = 0.05) # Define breaks for cosine scores
strata_COS <- setNames(
  lapply(seq_along(cos_breaks)[-length(cos_breaks)], function(i) c(cos_breaks[i], cos_breaks[i + 1])),
  paste0(sprintf("%.2f", cos_breaks[-length(cos_breaks)]), "-", sprintf("%.2f", cos_breaks[-1])))


## Function to stratify and count values, ignoring NAs
stratify_and_count <- function(df, column, strata) {
  counts <- sapply(strata, function(range) {
    sum(df[[column]] >= range[1] & df[[column]] < range[2], na.rm = TRUE) # Count values within each stratum
  })
  return(counts)
}

## Create lists to save individual df and strata counts
individual_dfs_COS <- list()
strata_counts_COS <- list()

## For each sample, create an individual df and strata counts
for (sample_name in unique(score_COS$ATTRIBUTE_Strain)) {
  individual_df <- score_COS %>%
    filter(ATTRIBUTE_Strain == sample_name) %>%
    select(id, Cosine_score)
  
  # Counting unique values in each stratum for Cosine score
  cosine_score_counts <- stratify_and_count(individual_df, "Cosine_score", strata_COS)
  
  # Create a dataframe for the counts
  counts_df <- data.frame(Stratum = names(cosine_score_counts), Count = cosine_score_counts)
  
  # Store count dataframe in the list
  strata_counts_COS[[sample_name]] <- counts_df
  
  # Assign the individual dataframe to the list
  individual_dfs_COS[[sample_name]] <- individual_df
}

## Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder
for(sample_name in names(individual_dfs_COS)) {
  write.csv(individual_dfs_COS[[sample_name]], file = paste0("Data/CSV_IQ/", sample_name, "_data_COS.csv"), row.names = FALSE)

  
  # Assign each dataframe to a variable in the R environment
  assign(paste0(sample_name, "_data_COS"), individual_dfs_COS[[sample_name]])
  assign(paste0(sample_name, "_strata_counts_COS"), strata_counts_COS[[sample_name]])
}

## name list (strains)
strains <- unique(Total_nodes$ATTRIBUTE_strain)

## df for summary
shannon_summary_cos <- data.frame(Strain = character(), Shannon_Index = numeric(), stringsAsFactors = FALSE)


## Iterate each strain
for (strain in strains) {
  # df name
  df_name <- paste0(strain, "_strata_counts_COS")
  if (exists(df_name)) {
    df_strata_counts <- get(df_name)
    
    # obtain total node for each strain from total_nodes
    total <- Total_nodes$Count[Total_nodes$ATTRIBUTE_strain == strain]

    if (total > 0) {
      if ("Count" %in% colnames(df_strata_counts)) {
        # Calculate the proportion (pi) for each stratum and add it as a new column
        df_strata_counts$pi <- df_strata_counts$Count / total
        
        # add a column with the total number
        df_strata_counts$Total <- total
        
        # calculate pi * ln(pi), handling the case of pi = 0
        df_strata_counts$pi_ln_pi <- ifelse(df_strata_counts$pi > 0, df_strata_counts$pi * log(df_strata_counts$pi), 0)
      } else {
        warning(paste("The column 'Count' does not exist in the dataframe", df_name))
        df_strata_counts$pi <- NA
        df_strata_counts$Total <- NA
        df_strata_counts$pi_ln_pi <- NA
      }
    } else {
      df_strata_counts$pi <- NA
      df_strata_counts$Total <- NA
      df_strata_counts$pi_ln_pi <- NA
    }
    # Calculate SI
    shannon_index <- -sum(df_strata_counts$pi_ln_pi, na.rm = TRUE)
    df_strata_counts$shannon_index <- shannon_index
    shannon_summary_cos <- rbind(shannon_summary_cos, data.frame(Strain = strain, Shannon_Index = shannon_index))
    df_strata_counts <- df_strata_counts[, c("Stratum", "Count", "pi", "Total", "pi_ln_pi", "shannon_index")]
    
  
    write.csv(df_strata_counts, file = paste0("Data/CSV_IQ/", df_name, "_SI_COS.csv"), row.names = FALSE)
  } else {
    warning(paste("El dataframe", df_name, "no existe."))
  }
}
write.csv(shannon_summary_cos, file = "Data/CSV_IQ/1_COS_shannon_index_summary.csv", row.names = FALSE)


# Annotated nodes per sample
Proportion_Annotated_nodes= merge(nodes_annotated_strain, Total_nodes, by.x = "ATTRIBUTE_strain", by.y = "ATTRIBUTE_strain")
Proportion_Annotated_nodes$Ans=Proportion_Annotated_nodes$Count.x/Proportion_Annotated_nodes$Count.y
Proportion_Annotated_nodes= Proportion_Annotated_nodes %>% select(-Count.x, -Count.y)

# Parameters union

IQ <- merge(Proportion_Annotated_nodes, shannon_summary_cos, by.x = "ATTRIBUTE_strain", by.y = "Strain")
IQ <- merge(IQ, shannon_summary_MQS, by.x = "ATTRIBUTE_strain", by.y = "Strain")
colnames(IQ)= c("ATTRIBUTE_Strain", "Ans", "Hcos", "Hmqs")

# Scaling
IQ$Hcos[is.na(IQ$Hcos)] <- 0
IQ$Hmqs[is.na(IQ$Hmqs)] <- 0
IQ$Ans[is.na(IQ$Ans)] <- 0
IQ$Hmqs <- (IQ$Hmqs) / (max(IQ$Hmqs))
IQ$Hcos <- (IQ$Hcos) / (max(IQ$Hcos))
IQ$IQ= (IQ$Ans+IQ$Hcos+IQ$Hmqs)/3
IQ = merge(IQ, Metadata, by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain")
IQ= IQ %>% select(IQ,Hmqs,Hcos, Ans, ATTRIBUTE_Bacteria, ATTRIBUTE_Strain)
IQ <- IQ %>%
  distinct()

num_unique_strains <- nrow(IQ)

# Create a custom cadet blue color gradient

# Create the plot with the new matte color palette and adjust the font size of the legend
IQ_plot = ggplot(IQ, aes(x = reorder(ATTRIBUTE_Strain, -IQ), y = IQ, fill = ATTRIBUTE_Bacteria)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "Black", width = 0.7, show.legend = TRUE) +  
  theme_classic() + 
  labs(x = "Strain",
       y = "IQ",
       fill = "Bacteria") +  # Legend label
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  
    axis.text.x = element_text(angle = 0, hjust = 1, size = 22),  
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 22),  
    legend.title = element_text(size = 22)  
  )+ scale_fill_viridis_d()

IQ_plot 


ggsave("Data/IQ_Plot/1_IQ_Plot.png", plot = IQ_plot, width = 15, height = 8)

# Distribution Plot
summary_data <- IQ %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(
    mean_IQ = mean(IQ, na.rm = TRUE),  # Calculate mean IQ
    sd_IQ = sd(IQ, na.rm = TRUE),      # Calculate standard deviation of IQ
    .groups = 'drop'                   # Drop grouping structure
  )

Distribution_plot = ggplot() +
  geom_point(data = IQ, aes(x = ATTRIBUTE_Bacteria, y = IQ), 
             shape = 1, color = "blue", size = 1.5, alpha = 0.6) +  # Puntos individuales
  geom_errorbar(data = summary_data, 
                aes(x = ATTRIBUTE_Bacteria, ymin = mean_IQ - sd_IQ, ymax = mean_IQ + sd_IQ), 
                width = 0.1, color = "black") +  # Barras de error para la media ± SD
  geom_point(data = summary_data, aes(x = ATTRIBUTE_Bacteria, y = mean_IQ), 
             shape = 1, color = "blue", size = 1) +  # Puntos de media
  labs(
    x = "", 
    y = "IQ") +  # Etiquetas de los ejes
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),    
    axis.text.y = element_text(size = 12),    
    axis.line = element_line(size = 1)        
  ) 

Distribution_plot

# Save the Distribution Plot as a PNG file
ggsave("Data/IQ_Plot/2_Distribution_plot.png", plot = Distribution_plot, width = 15, height = 8)

# Create Boxplot
Boxplot = ggplot(IQ, aes(x = ATTRIBUTE_Bacteria, y = IQ)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.1) +  # Boxplot without outliers
  geom_point(aes(color = "Individual Points"), shape = 1, size = 3, alpha = 0.6) +  # Individual points
  scale_color_manual(values = "blue") +  # Manual color for points
  labs(
    x = "Bacteria", 
    y = "IQ") +  # Axis labels
  theme_classic() +  # Minimalist theme
  theme(legend.position = "none")  # Remove legend

# Save the Boxplot as a PNG file
ggsave("Data/IQ_Plot/3_Boxplot.png", plot = Boxplot, width = 15, height = 8)

# Heatmap
heatmap = merge(Total_nodes_per_strain, Metadata, by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain")
heatmap = heatmap %>% select(id, ATTRIBUTE_Strain, ATTRIBUTE_Bacteria)

## Create a list of unique IDs for each combination of strain and bacteria
df_grouped <- heatmap %>% 
  group_by(ATTRIBUTE_Strain, ATTRIBUTE_Bacteria) %>% 
  summarise(id_list = list(unique(id)), .groups = 'drop')

## Create a nested list of IDs by bacteria and strain
id_lists <- split(df_grouped$id_list, df_grouped$ATTRIBUTE_Bacteria)

## Obtain the names of the bacteria and the strains
bacteria_names <- names(id_lists)

## Create a list of intersection matrices for each bacteria
intersect_matrices <- list()

for (bacteria in bacteria_names) {
  id_per_strain <- id_lists[[bacteria]]
  strain_names <- df_grouped$ATTRIBUTE_Strain[df_grouped$ATTRIBUTE_Bacteria == bacteria]
  n_strains <- length(strain_names)
  
  # Intersection matrix for the specific bacteria
  intersect_matrix <- matrix(0, nrow = n_strains, ncol = n_strains, dimnames = list(strain_names, strain_names))
  
  for (i in 1:n_strains) {
    for (j in 1:n_strains) {
      intersect_matrix[i, j] <- length(intersect(id_per_strain[[i]], id_per_strain[[j]]))
    }
  }
  
  # Store the matrix in a list
  intersect_matrices[[bacteria]] <- intersect_matrix
}

## Convert matrices to long format
intersect_matrix_melt <- do.call(rbind, lapply(names(intersect_matrices), function(bacteria) {
  mat_melt <- as.data.frame(as.table(intersect_matrices[[bacteria]]))
  mat_melt$Bacteria <- bacteria
  return(mat_melt)
}))

## Plot the heatmap with facets by bacteria
Heatmap <- ggplot(intersect_matrix_melt, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("white", "lightblue", "blue", "darkblue")) +
  labs(x = "", y = "", fill = "Shared nodes number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(size = 4)) 

Heatmap
## Save the Distribution Plot as a PNG file
ggsave("Data/IQ_Plot/4_Heatmap.png", plot = Heatmap, width = 15, height = 8)


# Network Analysis
library(igraph)
library(ggraph)
library(proxy)
library(visNetwork)

## Convert to matrix
data_matrix <- as.matrix(IQ[, c("Hmqs", "Hcos", "Ans")])

## Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

## Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ$ATTRIBUTE_Strain

## Define threshold
threshold <- 0.995

## Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

## Convert to numeric matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

## Create the graph from the adjacency matrix
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

## Add ATTRIBUTE_Bacteria as a vertex attribute
V(g)$bacteria <- IQ$ATTRIBUTE_Bacteria

## Export the graph as a GraphML file
write_graph(g, file = "Data/IQ_Plot/Similarity.graphml", format = "graphml")

## Create a color vector based on ATTRIBUTE_Bacteria
IQ$ATTRIBUTE_Bacteria <- as.factor(IQ$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ$ATTRIBUTE_Bacteria)))[as.numeric(IQ$ATTRIBUTE_Bacteria)]

## Visualize the network
plot(g, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 5, 
     vertex.color = color_vector,  # Use the color vector
     main = "")

# Summary of top 10 lowest IQ
IQ_Summary_top <- IQ %>%
  arrange(IQ) %>%
  slice_head(n = 10)

## IQ top plot

IQ_top_plot = ggplot(IQ_Summary_top, aes(x = reorder(ATTRIBUTE_Strain, -IQ), y = IQ, fill = ATTRIBUTE_Bacteria)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "Black", width = 0.7, show.legend = TRUE) +  
  theme_classic() + 
  labs(x = "",
       y = "IQ",
       fill = "Bacteria") +  # Legend label
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 14),  
    axis.title = element_text(face = "bold"),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 14)  
  )+ scale_fill_viridis_d()

IQ_top_plot 


IQ_top_plot
## Save the plot
ggsave("Data/IQ_Plot/5_IQ_top_Plot.png", plot = IQ_top_plot, width = 15, height = 8)

## Distribution top plot data
summary_top_data <- IQ_Summary_top %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(
    mean_IQ = mean(IQ, na.rm = TRUE),
    sd_IQ = sd(IQ, na.rm = TRUE),
    .groups = 'drop'
  )

## Distribution top plot
Distribution_top_plot = ggplot() +
  geom_point(data = IQ_Summary_top, aes(x = ATTRIBUTE_Bacteria, y = IQ), 
             shape = 1, color = "blue", size = 1.5, alpha = 0.6) +  # Points for individual IQ values
  geom_errorbar(data = summary_top_data, 
                aes(x = ATTRIBUTE_Bacteria, ymin = mean_IQ - sd_IQ, ymax = mean_IQ + sd_IQ), 
                width = 0.1, color = "Black") +  # Error bars for mean IQ
  geom_point(data = summary_top_data, aes(x = ATTRIBUTE_Bacteria, y = mean_IQ), 
             shape = 1, color = "blue", size = 1) +  # Points for mean IQ
  labs(
    x = "Bacteria", 
    y = "IQ") +
  theme_classic() 

Distribution_top_plot
# Save the distribution plot
ggsave("Data/IQ_Plot/6_Distribution_top_plot.png", plot = Distribution_top_plot, width = 15, height = 8)

# Heatmap top

heatmap_top = merge(heatmap, IQ_Summary_top, by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain")
heatmap_top = heatmap_top %>% select(id, ATTRIBUTE_Strain, ATTRIBUTE_Bacteria.x)
colnames(heatmap_top) = c("id", "ATTRIBUTE_strain", "ATTRIBUTE_Bacteria")

## Create a list of unique IDs for each strain and bacteria combination
df_grouped_top <- heatmap_top %>% 
  group_by(ATTRIBUTE_strain, ATTRIBUTE_Bacteria) %>% 
  summarise(id_list = list(unique(id)))

## Create a nested list of IDs by bacteria and strain
id_lists_top <- split(df_grouped_top$id_list, df_grouped_top$ATTRIBUTE_Bacteria)

## Get the names of the bacteria and strains
bacteria_names_top <- names(id_lists_top)

## Create a list of intersection matrices for each bacteria
intersect_matrices_top <- list()

for (bacteria in bacteria_names_top) {
  id_per_strain_top <- id_lists_top[[bacteria]]
  strain_names_top <- df_grouped_top$ATTRIBUTE_strain[df_grouped_top$ATTRIBUTE_Bacteria == bacteria]
  n_strains_top <- length(strain_names_top)
  
  # Intersection matrix for the specific bacteria
  intersect_matrix <- matrix(0, nrow = n_strains_top, ncol = n_strains_top, dimnames = list(strain_names_top, strain_names_top))
  for (i in 1:n_strains_top) {
    for (j in 1:n_strains_top) {
      intersect_matrix[i, j] <- length(intersect(id_per_strain_top[[i]], id_per_strain_top[[j]]))
    }
  }
  
  # Store the matrix in a list
  intersect_matrices_top[[bacteria]] <- intersect_matrix
}

## Convert the matrices to long format
intersect_matrix_melt_top <- do.call(rbind, lapply(names(intersect_matrices_top), function(bacteria) {
  mat_melt <- as.data.frame(as.table(intersect_matrices_top[[bacteria]]))  # Updated to intersect_matrices_top
  mat_melt$Bacteria <- bacteria
  return(mat_melt)
}))

## Plot the heatmap with facets by bacteria
Heatmap_top <- ggplot(intersect_matrix_melt_top, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("white", "lightblue", "blue", "darkblue")) +
  labs(x = "", y = "", fill = "Shared nodes number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) 
Heatmap_top

# Save the distribution plot
ggsave("Data/IQ_Plot/7_Heatmap_top_plot.png", plot = Heatmap_top, width = 15, height = 8)


# Red analysis top 10

# Convert to matrix
data_matrix <- as.matrix(IQ_Summary_top[, c("Hmqs", "Hcos", "Ans")])

# Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

# Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ_Summary_top$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ_Summary_top$ATTRIBUTE_Strain

# Define threshold
threshold <- 0.995

# Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

# Convert to numerical matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

# Create the graph from the adjacency matrix
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Ensure ATTRIBUTE_Bacteria is character or factor
IQ_Summary_top$ATTRIBUTE_Bacteria <- as.character(IQ_Summary_top$ATTRIBUTE_Bacteria)

# Add ATTRIBUTE_Bacteria as a vertex attribute
V(g)$bacteria <- IQ_Summary_top$ATTRIBUTE_Bacteria

# Export the graph as a GraphML file
write_graph(g, file = "Data/IQ_Plot/Top_graphlm.graphml", format = "graphml")

# Create a color vector based on ATTRIBUTE_Bacteria
IQ_Summary_top$ATTRIBUTE_Bacteria <- as.factor(IQ_Summary_top$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ_Summary_top$ATTRIBUTE_Bacteria)))[as.numeric(IQ$ATTRIBUTE_Bacteria)]

# Visualize the network
plot(g, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 30, 
     vertex.color = color_vector,  # Use the color vector
     main = "")

# Proportions plot

# Count how many ATTRIBUTE_strain there are for each ATTRIBUTE_Bacteria in the original dataframe
strain_count_per_bacteria <- Metadata %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(n_strains = n_distinct(ATTRIBUTE_Strain))

# Count how many ATTRIBUTE_strain there are for each ATTRIBUTE_Bacteria in the top dataframe
strain_count_per_bacteria_top <- IQ_Summary_top %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(n_strains_top = n_distinct(ATTRIBUTE_Strain))

# Join both dataframes to calculate the proportion
proportion_counts <- strain_count_per_bacteria %>%
  left_join(strain_count_per_bacteria_top, by = "ATTRIBUTE_Bacteria") %>%
  mutate(n_strains_top = replace_na(n_strains_top, 0),  # Replace NA with 0
         proportion = n_strains_top / n_strains * 100)
proportion_counts$n_strain_top_proportion<- proportion_counts$n_strains_top/10*100

# Plot the proportions
PPlot<- ggplot(proportion_counts, aes(x = reorder(ATTRIBUTE_Bacteria, -proportion), y = proportion, fill = proportion)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "cadetblue") + 
  labs(
    x = "",
    y = "Strain (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8))

ggsave("Data/IQ_Plot/8_Proportion_plot.png", plot = PPlot, width = 15, height = 8)


Bb_plot <- ggplot(proportion_counts, aes(x = proportion, y = n_strains_top, size = n_strains, color = ATTRIBUTE_Bacteria, label = ATTRIBUTE_Bacteria)) +
  geom_point(alpha = 0.7) + 
  geom_text(vjust = -2, hjust = 0.5, check_overlap = TRUE, size = 3) +  
  scale_size_continuous(range = c(3, 10)) +  
  scale_color_brewer(palette = "Set2") +  
  labs(x = "Relative proportion (%) [Number of strain in top / total strain]", 
       y = "Proportion in top (%) [Number of strain in top / top (n)]", 
       size = "Total strains") +  # Eliminamos el título de la leyenda del color
  theme_classic() +  
  theme(legend.position = "right",  
        axis.title = element_text(size = 12),  
        axis.text = element_text(size = 10),  
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  guides(color = "none")  # Elimina la leyenda de color

ggsave("Data/IQ_Plot/9_Proportion_Bubble_plot.png", plot = Bb_plot, width = 15, height = 8)

# Bini integration

# BiNi result

BiNI_results <- read_excel("Data/BiNI_results.xlsx")


IQ_BiNI <- merge(IQ, BiNI_results, by.x = "ATTRIBUTE_Strain", by.y = "Strain")
IQ_BiNI_long <- IQ_BiNI %>% select(ATTRIBUTE_Bacteria, ATTRIBUTE_Strain, IQ, BiNI)

# IQ-BiNi plot
IQ_BiNI_long <- IQ_BiNI_long %>%
  pivot_longer(
    cols = -c(ATTRIBUTE_Strain, ATTRIBUTE_Bacteria),
    names_to = "Parameter",    
    values_to = "Value"        
  )

# Create a custom cadet blue color gradient
palette <- colorRampPalette(c("#5F9EA0", "#4682B4", "#2F4F4F"))(length(unique(IQ_BiNI_long$ATTRIBUTE_Bacteria)))

# Create the plot with dynamically generated aquamarine colors
## Create the plot
IQ_BiNI_Plot <- ggplot(IQ_BiNI_long, aes(x = ATTRIBUTE_Strain, y = Value, color = ATTRIBUTE_Bacteria, shape = Parameter)) +
  geom_point(size = 3, alpha = 0.8) +  # Add transparency for better overlap handling
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_classic(base_size = 14) +  # Base font size increased for better readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels for better fit
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),  # Add margin to the x label
    axis.title.y = element_text(size = 12, margin = margin(r = 10))
  ) +
  ylim(0, 1) +
  scale_shape_manual(values = c(16, 17)) +  # Shape 16 for IQ, shape 17 for BiNI
  scale_color_manual(values = palette) +  # Use the defined aquamarine palette
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
  scale_color_viridis_d() 

IQ_BiNI_Plot

ggsave("Data/IQ_Plot/10_IQ_BiNI_plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)


## Create the scatter plot of BiNI vs IQ for each strain
IQ_BiNI_Plot <- ggplot(IQ_BiNI, aes(x = BiNI, y = IQ, color = ATTRIBUTE_Bacteria, shape = ATTRIBUTE_Strain)) +
  geom_point(size = 3) +
  labs(title = "",
       x = "BiNI",
       y = "IQ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylim(0, 1) +
  xlim(0, 1) +  # Ensure both axes are on the same scale
  scale_shape_manual(values = 1:length(unique(IQ_BiNI$ATTRIBUTE_Strain))) +  
  scale_fill_viridis_d()  

IQ_BiNI_Plot

ggsave("Data/IQ_Plot/11_IQ_BiNI_s_Plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)


# scatter plot IQ_BiNI
# Define a color palette

palette <- brewer.pal(n = length(unique(IQ_BiNI$ATTRIBUTE_Bacteria)), name = "Set2")
correlation <- cor(IQ_BiNI$BiNI, IQ_BiNI$IQ, use = "complete.obs")  


IQ_BiNI_SPlot <- ggplot(IQ_BiNI, aes(x = IQ, y = BiNI, color = ATTRIBUTE_Bacteria)) +
  geom_point(size = 4, alpha = 0.7) +  # Adjusting point transparency and size
  geom_smooth(method = "lm", linetype = "dashed", color = "black", fill = "lightgray") +  # Confidence interval
  labs(title = "",
       x = "IQ",
       y = "BiNI") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),  # Increase axis text size
    axis.title = element_text(size = 14),  # Increase axis title size
    plot.title = element_text(hjust = 0.5, size = 16)  # Center and increase plot title size
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("R: ", round(correlation, 2)), 
           hjust = 3, vjust = 10, size = 5, color = "black") +
  scale_color_viridis_d() 
IQ_BiNI_SPlot


ggsave("Data/IQ_Plot/11_IQ_BiNI_SPlot.png", plot = IQ_BiNI_SPlot, width = 15, height = 8)


# Convert to matrix
IQ_BiNI_f <- IQ_BiNI %>% select(-BiNI)
data_matrix <- as.matrix(IQ_BiNI_f[, c("Hmqs", "Hcos", "Ans")])

# Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

# Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ_BiNI_f$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ_BiNI_f$ATTRIBUTE_Strain

# Define threshold
threshold <- 0.995

# Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

# Convert to numerical matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

# Create the graph from the adjacency matrix
IQ_BiNI_G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Ensure ATTRIBUTE_Bacteria is character or factor
IQ_BiNI_f$ATTRIBUTE_Bacteria <- as.character(IQ_BiNI_f$ATTRIBUTE_Bacteria)

# Add ATTRIBUTE_Bacteria as a vertex attribute
V(IQ_BiNI_G)$bacteria <- IQ_BiNI_f$ATTRIBUTE_Bacteria

# Export the graph as a GraphML file
write_graph(IQ_BiNI_G, file = "Data/IQ_Plot/IQ_BiNi_graphml.graphml", format = "graphml")

# Create a color vector based on ATTRIBUTE_Bacteria
IQ_BiNI_f$ATTRIBUTE_Bacteria <- as.factor(IQ_BiNI_f$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ_BiNI_f$ATTRIBUTE_Bacteria)))[as.numeric(IQ_BiNI_f$ATTRIBUTE_Bacteria)]

## Visualize the network
plot(IQ_BiNI_G, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 30, 
     vertex.color = color_vector,  # Use the color vector
     main = "")







# NPAtlas integration

# Count total nodes per strain
Total_nodes <- table(Total_nodes_per_strain$ATTRIBUTE_Strain)
Total_nodes <- as.data.frame(Total_nodes)
colnames(Total_nodes) <- c("ATTRIBUTE_strain", "Count")


NPAtlas_data <- read_excel("NPAtlas_data.xlsx")
NPAtlas_data= NPAtlas_data %>% select(compound_name, compound_inchi, compound_smiles, origin_type, genus, origin_species)

# Extract index inputs
score_MQS <- df_nodes %>%
  select(id, MQScore, ATTRIBUTE_Strain, Smiles, INCHI, Compound_Name) %>%
  filter(!is.na(MQScore) & MQScore != "") # Extract MSQ and remove NAs

# canonizise smiles

NPAtlas_data$Canonical_Smiles <- sapply(NPAtlas_data$compound_smiles, function(smiles) {
  mol <- parse.smiles(smiles)[[1]]
  get.smiles(mol)
})
NPAtlas_data= NPAtlas_data %>% select(-compound_smiles)

score_MQS$Canonical_Smiles= sapply(score_MQS$Smiles, function(smiles){
  mol <- parse.smiles(smiles)
  if (is.null(mol[[1]])) {
    return(NA)  
  } else {
    return(get.smiles(mol[[1]]))
  }
})

score_MQS= score_MQS %>% select(-Smiles)

# get NPAtlas anotations

#By smiles
MQS_score_NPAtlas_smiles <- left_join(score_MQS, NPAtlas_data, by = "Canonical_Smiles", 
                                      relationship = "many-to-many")%>% filter(!is.na(origin_type))
# By inchi
MQS_score_NPAtlas_inchi <- left_join(score_MQS, NPAtlas_data, 
                                     by = c("INCHI" = "compound_inchi"), 
                                     relationship = "many-to-many") %>%
  filter(!is.na(origin_type))

#By name
MQS_score_NPAtlas_name <- left_join(score_MQS, NPAtlas_data, 
                                    by = c("Compound_Name" = "compound_name"), 
                                    relationship = "many-to-many") %>%
  filter(!is.na(origin_type))

# COMBINE
MQS_score_NPAtlas <- bind_rows(MQS_score_NPAtlas_smiles, MQS_score_NPAtlas_inchi, MQS_score_NPAtlas_name) %>% 
  distinct(id, .keep_all = TRUE) %>%
  select(id, MQScore, ATTRIBUTE_Strain)  %>% 
  separate_rows(ATTRIBUTE_Strain, sep = "\\,") %>%  
  group_by(ATTRIBUTE_Strain)

MQS_score_NPAtlas<-merge(MQS_score_NPAtlas, Metadata, 
                         by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain", all.x = TRUE)

MQS_score_NPAtlas <- MQS_score_NPAtlas %>%
  select(id, MQScore, ATTRIBUTE_Strain) %>%
  filter(!is.na(ATTRIBUTE_Strain))

# Count nodes with annotated strain information
nodes_annotated_strain_NPAtlas <- table(MQS_score_NPAtlas$ATTRIBUTE_Strain)
nodes_annotated_strain_NPAtlas <- as.data.frame(nodes_annotated_strain_NPAtlas)
colnames(nodes_annotated_strain_NPAtlas) <- c("ATTRIBUTE_strain", "Count")

## Define strata for MQScore
MQS_breaks <- seq(0.00, 1, by = 0.05)
strata_MQS <- setNames(
  lapply(seq_along(MQS_breaks)[-length(MQS_breaks)], function(i) c(MQS_breaks[i], MQS_breaks[i + 1])),
  paste0(sprintf("%.2f", MQS_breaks[-length(MQS_breaks)]), "-", sprintf("%.2f", MQS_breaks[-1])))

## Function to stratify and count values, ignoring NAs
stratify_and_count_unique <- function(df, column, strata) {
  unique_df <- df %>% 
    group_by(id) %>% 
    slice(1) %>%  # Count only the first occurrence of each unique identifier (id) in the specified column
    ungroup()
  
  # Count unique values for each stratum
  counts <- sapply(strata, function(range) {
    sum(unique_df[[column]] >= range[1] & unique_df[[column]] < range[2], na.rm = TRUE)
  })
  
  return(counts)
}


# Initialize empty lists to store results
strata_counts <- list()
individual_dfs <- list()

# Loop over each unique sample
for (sample_name in unique(MQS_score_NPAtlas$ATTRIBUTE_Strain)) {
  # Create a temporary dataframe for the individual sample
  individual_df <- MQS_score_NPAtlas %>%
    filter(ATTRIBUTE_Strain == sample_name) %>%
    select(id, MQScore)
  
  # Count unique values in each stratum for MQScore
  MQScore_counts <- stratify_and_count_unique(individual_df, "MQScore", strata_MQS)
  
  # Create a dataframe with consistent strata names, ensuring that all strata are represented
  counts_df <- data.frame(
    Stratum = names(strata_MQS),
    MQScore_Count = rep(0, length(strata_MQS))  # Initialize with zeros
  )
  
  # Populate the counts for the strata that exist
  for (i in seq_along(names(strata_MQS))) {
    stratum_name <- names(strata_MQS)[i]
    if (stratum_name %in% names(MQScore_counts)) {
      counts_df$MQScore_Count[i] <- MQScore_counts[[stratum_name]]
    }
  }
  
  # Store the results in lists
  strata_counts[[sample_name]] <- counts_df
  individual_dfs[[sample_name]] <- individual_df
}


## Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder
for(sample_name in names(individual_dfs)) {
  # Save individual dataframe to CSV
  write.csv(individual_dfs[[sample_name]], 
            file = paste0("Data/CSV_IQNPA/", sample_name, "_NPAtlas_data_MQS.csv"), 
            row.names = FALSE)
  
  # Save strata counts to CSV
  write.csv(strata_counts[[sample_name]], 
            file = paste0("Data/CSV_IQNPA/", sample_name, "_NPAtlas_strata_counts.csv"), 
            row.names = FALSE)
  
  # Assign each dataframe to a variable in the R environment
  assign(paste0(sample_name, "_NPAtlas_data"), 
         individual_dfs[[sample_name]], 
         envir = .GlobalEnv)  # Specify environment to ensure it's assigned globally
  assign(paste0(sample_name, "_NPAtlas_strata_counts"), 
         strata_counts[[sample_name]], 
         envir = .GlobalEnv)
}

## Create a data frame to store strain names
strains <- data.frame(Total_nodes$ATTRIBUTE_strain)
colnames(strains) <- c("strain")
strains <- strains$strain

## Create an empty dataframe for the summary
shannon_summary_MQS <- data.frame(Strain = character(), Shannon_Index = numeric(), stringsAsFactors = FALSE)

## Iterate over each strain
for (strain in strains) {
  # Name of the corresponding dataframe
  df_name <- paste0(strain, "_NPAtlas_strata_counts")
  
  # Check if the dataframe exists
  if (exists(df_name)) {
    df_strata_counts <- get(df_name)
    
    # Obtain the total count for this strain from Total_nodes
    total <- Total_nodes$Count[Total_nodes$ATTRIBUTE_strain == strain]
    
    # Check that the total is greater than zero
    if (total > 0) {
      # Calculate the proportion (pi) for each stratum and add it as a new column
      df_strata_counts$pi <- df_strata_counts$MQScore_Count / total
      
      # Add a column with the total count
      df_strata_counts$Total <- total
      
      # Calculate pi * ln(pi), handling the case where pi = 0
      df_strata_counts$pi_ln_pi <- ifelse(df_strata_counts$pi > 0, df_strata_counts$pi * log(df_strata_counts$pi), 0)
    } else {
      # If the total is zero, assign NA to the new columns
      df_strata_counts$pi <- NA
      df_strata_counts$Total <- NA
      df_strata_counts$pi_ln_pi <- NA
    }
    
    # Calculate the Shannon index as -sum(pi * ln(pi))
    shannon_index <- -sum(df_strata_counts$pi_ln_pi, na.rm = TRUE)
    
    # Append the Shannon index to the summary
    shannon_summary_MQS <- rbind(shannon_summary_MQS, data.frame(Strain = strain, Shannon_Index = shannon_index))
    
    # Clean column names to remove any whitespace
    colnames(df_strata_counts) <- trimws(colnames(df_strata_counts))
    
    # Check if the dataframe has rows
    if (nrow(df_strata_counts) > 0) {
      # Define required columns
      required_columns <- c("Stratum", "MQScore_Count", "pi", "Total", "pi_ln_pi", "shannon_index")
      
      # Check for missing columns
      missing_columns <- setdiff(required_columns, colnames(df_strata_counts))
      if (length(missing_columns) > 0) {
        warning(paste("Faltan las siguientes columnas en el dataframe:", paste(missing_columns, collapse = ", ")))
      } else {
        # Select only the relevant columns
        df_strata_counts <- df_strata_counts[, required_columns]
        
        # Add the Shannon index as a new column (same for all strata)
        df_strata_counts$shannon_index <- shannon_index
        
        # Save the updated dataframe to a CSV file
        write.csv(df_strata_counts, file = paste0("Data/CSV_IQ/", df_name, "_SI_MQS.csv"), row.names = FALSE)
      }
    } else {
      warning("El dataframe df_NPAtlas_strata_counts está vacío.")
    }
  } else {
    warning(paste("The dataframe", df_name, "does not exist."))
  }
}
warnings()

## Save the overall Shannon summary to a CSV file
write.csv(shannon_summary_MQS, file = "Data/CSV_IQNPA/1_Shannon_NPAtlas_Summary_MQS.csv", row.names = FALSE)


# Extract index input: COS_score

# Get cosine scores where node2 matches the id in score_MQS
score_COS_node1 <- df_edges %>%
  select(node1, cosine_score, node2) %>%
  semi_join(MQS_score_NPAtlas, by = c("node2" = "id")) # Retrieve cosine values from node to node

# Get complementary cosine scores where node1 matches the id in score_MQS
score_COS_node2 <- df_edges %>%
  select(node2, cosine_score, node1) %>%
  semi_join(MQS_score_NPAtlas, by = c("node1" = "id")) # Retrieve complementary cosine values

## Combine interactions into a single dataframe
combined_interactions <- bind_rows(
  score_COS_node1 %>% mutate(direction = "Score_COS_node1"),
  score_COS_node2 %>% mutate(direction = "Score_COS_node2")
) # Add a new column called direction to indicate the source node of the interaction.

# Create a combined node column based on direction
combined_interactions <- combined_interactions %>% mutate(
  combined_nodes = ifelse(direction == "Score_COS_node1",
                          paste(node1, node2, sep = "/"),
                          paste(node2, node1, sep = "/"))
) # If “Score_COS_node1”, concatenate node1 and node2; else, concatenate node2 and node1.

# Separate combined nodes into nodeA and nodeB
combined_interactions <- combined_interactions %>%
  separate(combined_nodes, into = c("nodeA", "nodeB"), sep = "/") # Assign nodeA as the node of interest.

## Remove rows where nodeA and nodeB are the same (singlets)
combined_interactions <- combined_interactions[combined_interactions$nodeA != combined_interactions$nodeB, ]

# Select relevant columns for cosine scores
score_COS <- combined_interactions %>%
  select(nodeA, cosine_score) # Filter cosine values for the node of interest.

# Merge cosine scores with node metadata
score_COS <- merge(score_COS, df_nodes, by.x = "nodeA", by.y = "id", all = TRUE)  
# Note: MQS values are repeated when there are multiple COS scores, as intended.

# Clean up score_COS dataframe
score_COS <- score_COS %>% select(nodeA, cosine_score, ATTRIBUTE_Strain) %>% 
  filter(!is.na(cosine_score) & cosine_score != "") %>% 
  separate_rows(ATTRIBUTE_Strain, sep = "\\,") %>%  
  group_by(ATTRIBUTE_Strain)

# Rename columns for clarity
colnames(score_COS) <- c("id", "Cosine_score", "ATTRIBUTE_Strain")

# Merge with metadata to include additional information
score_COS <- merge(score_COS, Metadata, 
                   by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain", all.x = TRUE)

# Final selection of relevant columns
score_COS <- score_COS %>% select(id, Cosine_score, ATTRIBUTE_Strain)

## Define strata for cosine score
cos_breaks <- seq(0.70, 1.00, by = 0.05) # Define breaks for cosine scores
strata_COS <- setNames(
  lapply(seq_along(cos_breaks)[-length(cos_breaks)], function(i) c(cos_breaks[i], cos_breaks[i + 1])),
  paste0(sprintf("%.2f", cos_breaks[-length(cos_breaks)]), "-", sprintf("%.2f", cos_breaks[-1])))


## Function to stratify and count values, ignoring NAs
stratify_and_count <- function(df, column, strata) {
  counts <- sapply(strata, function(range) {
    sum(df[[column]] >= range[1] & df[[column]] < range[2], na.rm = TRUE) # Count values within each stratum
  })
  return(counts)
}

## Create lists to save individual df and strata counts
individual_dfs_COS <- list()
strata_counts_COS <- list()

## For each sample, create an individual df and strata counts
for (sample_name in unique(score_COS$ATTRIBUTE_Strain)) {
  individual_df <- score_COS %>%
    filter(ATTRIBUTE_Strain == sample_name) %>%
    select(id, Cosine_score)
  
  # Counting unique values in each stratum for Cosine score
  cosine_score_counts <- stratify_and_count(individual_df, "Cosine_score", strata_COS)
  
  # Create a dataframe for the counts
  counts_df <- data.frame(Stratum = names(cosine_score_counts), Count = cosine_score_counts)
  
  # Store count dataframe in the list
  strata_counts_COS[[sample_name]] <- counts_df
  
  # Assign the individual dataframe to the list
  individual_dfs_COS[[sample_name]] <- individual_df
}

## Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder
for(sample_name in names(individual_dfs_COS)) {
  write.csv(individual_dfs_COS[[sample_name]], file = paste0("Data/CSV_IQNPA/", sample_name, "_data_COS.csv"), row.names = FALSE)
  
  
  # Assign each dataframe to a variable in the R environment
  assign(paste0(sample_name, "_data_COS"), individual_dfs_COS[[sample_name]])
  assign(paste0(sample_name, "_strata_counts_COS"), strata_counts_COS[[sample_name]])
}

## name list (strains)
strains <- unique(Total_nodes$ATTRIBUTE_strain)

## df for summary
shannon_summary_cos <- data.frame(Strain = character(), Shannon_Index = numeric(), stringsAsFactors = FALSE)

# Iterate through each strain
for (strain in strains) {
  # Create dataframe name
  df_name <- paste0(strain, "_strata_counts_COS")
  
  if (exists(df_name)) {
    # Retrieve the dataframe
    df_strata_counts <- get(df_name)
    
    # Obtain total node for each strain from total_nodes
    total <- Total_nodes$Count[Total_nodes$ATTRIBUTE_strain == strain]
    
    if (total > 0) {
      if ("Count" %in% colnames(df_strata_counts)) {
        # Calculate the proportion (pi) for each stratum and add it as a new column
        df_strata_counts$pi <- df_strata_counts$Count / total
        
        # Add a column with the total number
        df_strata_counts$Total <- total
        
        # Calculate pi * ln(pi), handling the case of pi = 0
        df_strata_counts$pi_ln_pi <- ifelse(df_strata_counts$pi > 0, df_strata_counts$pi * log(df_strata_counts$pi), 0)
      } else {
        warning(paste("The column 'Count' does not exist in the dataframe", df_name))
        df_strata_counts$pi <- NA
        df_strata_counts$Total <- NA
        df_strata_counts$pi_ln_pi <- NA
      }
    } else {
      df_strata_counts$pi <- NA
      df_strata_counts$Total <- NA
      df_strata_counts$pi_ln_pi <- NA
    }
    
    # Calculate Shannon Index (SI)
    shannon_index <- -sum(df_strata_counts$pi_ln_pi, na.rm = TRUE)
    df_strata_counts$shannon_index <- shannon_index
    
    # Append the strain and Shannon index to the summary
    shannon_summary_cos <- rbind(shannon_summary_cos, data.frame(Strain = strain, Shannon_Index = shannon_index))
    
    # Reorder columns and save the dataframe
    df_strata_counts <- df_strata_counts[, c("Stratum", "Count", "pi", "Total", "pi_ln_pi", "shannon_index")]
    write.csv(df_strata_counts, file = paste0("Data/CSV_IQNPA/", df_name, "_SI_COS.csv"), row.names = FALSE)
    
  } else {
    # If dataframe does not exist, create a dataframe with Count = 0 and Shannon index = 0
    df_strata_counts <- data.frame(Stratum = NA, Count = 0, pi = NA, Total = 0, pi_ln_pi = NA, shannon_index = 0)
    
    # Append Shannon index = 0 to the summary
    shannon_summary_cos <- rbind(shannon_summary_cos, data.frame(Strain = strain, Shannon_Index = 0))
    
    # Save the empty dataframe
    write.csv(df_strata_counts, file = paste0("Data/CSV_IQNPA/", df_name, "_SI_COS.csv"), row.names = FALSE)
    
    # Display a warning about the missing dataframe
    warning(paste("El dataframe", df_name, "no existe. Se generó un dataframe con conteo 0 e índice Shannon 0."))
  }
} 

# Save the Shannon summary
write.csv(shannon_summary_cos, file = "Data/CSV_IQNPA/1_COS_shannon_index_summary.csv", row.names = FALSE)
shannon_summary_cos <- shannon_summary_cos %>%
  distinct()

# Annotated nodes per sample
Proportion_Annotated_nodes= merge(nodes_annotated_strain_NPAtlas, Total_nodes, by.x = "ATTRIBUTE_strain", by.y = "ATTRIBUTE_strain")
Proportion_Annotated_nodes$Ans=Proportion_Annotated_nodes$Count.x/Proportion_Annotated_nodes$Count.y
Proportion_Annotated_nodes= Proportion_Annotated_nodes %>% select(-Count.x, -Count.y)

# Parameters union

IQ <- merge(Proportion_Annotated_nodes, shannon_summary_cos, by.x = "ATTRIBUTE_strain", by.y = "Strain", all = TRUE)
IQ <- merge(IQ, shannon_summary_MQS, by.x = "ATTRIBUTE_strain", by.y = "Strain", all = TRUE)
colnames(IQ)= c("ATTRIBUTE_Strain", "Ans", "Hcos", "Hmqs")

# Scaling
IQ$Hcos[is.na(IQ$Hcos)] <- 0
IQ$Hmqs[is.na(IQ$Hmqs)] <- 0
IQ$Ans[is.na(IQ$Ans)] <- 0
IQ$Hmqs <- (IQ$Hmqs) / (max(IQ$Hmqs))
IQ$Hcos <- (IQ$Hcos) / (max(IQ$Hcos))
IQ$IQ= (IQ$Ans+IQ$Hcos+IQ$Hmqs)/3
IQ = merge(IQ, Metadata, by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain")
IQ= IQ %>% select(IQ,Hmqs,Hcos, Ans, ATTRIBUTE_Bacteria, ATTRIBUTE_Strain)
IQ_NPA <- IQ %>%
  distinct()


# IQ plot
IQ_plot = ggplot(IQ_NPA, aes(x = reorder(ATTRIBUTE_Strain, -IQ), y = IQ, fill = ATTRIBUTE_Bacteria)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "Black", width = 0.7, show.legend = TRUE) +  
  theme_classic() + 
  labs(x = "Strain",
       y = "IQNP",
       fill = "Bacteria") +  # Legend label
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  
    axis.text.x = element_text(angle = 0, hjust = 1, size = 22),  
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 22),  
    legend.title = element_text(size = 22)  
  )+ scale_fill_viridis_d()


IQ_plot
ggsave("Data/IQNPA_Plot/1_IQ_Plot.png", plot = IQ_plot, width = 15, height = 8)

# Distribution Plot
summary_data <- IQ_NPA %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(
    mean_IQ = mean(IQ, na.rm = TRUE),  # Calculate mean IQ
    sd_IQ = sd(IQ, na.rm = TRUE),      # Calculate standard deviation of IQ
    .groups = 'drop'                   # Drop grouping structure
  )

Distribution_plot = ggplot() +
  geom_point(data = IQ, aes(x = ATTRIBUTE_Bacteria, y = IQ), 
             shape = 1, color = "blue", size = 1.5, alpha = 0.6) +  # Puntos individuales
  geom_errorbar(data = summary_data, 
                aes(x = ATTRIBUTE_Bacteria, ymin = mean_IQ - sd_IQ, ymax = mean_IQ + sd_IQ), 
                width = 0.1, color = "black") +  # Barras de error para la media ± SD
  geom_point(data = summary_data, aes(x = ATTRIBUTE_Bacteria, y = mean_IQ), 
             shape = 1, color = "blue", size = 1) +  # Puntos de media
  labs(
    x = "", 
    y = "IQ") +  # Etiquetas de los ejes
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),    
    axis.text.y = element_text(size = 12),    
    axis.line = element_line(size = 1)        
  ) 
Distribution_plot
# Save the Distribution Plot as a PNG file
ggsave("Data/IQNPA_Plot/2_Distribution_plot.png", plot = Distribution_plot, width = 15, height = 8)

# Create Boxplot
Boxplot = ggplot(IQ_NPA, aes(x = ATTRIBUTE_Bacteria, y = IQ)) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.1) +  # Boxplot without outliers
  geom_point(aes(color = "Individual Points"), shape = 1, size = 3, alpha = 0.6) +  # Individual points
  scale_color_manual(values = "blue") +  # Manual color for points
  labs(
    x = "", 
    y = "IQ") +  # Axis labels
  theme_classic() +  # Minimalist theme
  theme(legend.position = "none")  # Remove legend

Boxplot
# Save the Boxplot as a PNG file
ggsave("Data/IQNPA_Plot/3_Boxplot.png", plot = Boxplot, width = 15, height = 8)

# Heatmap
heatmap = merge(Total_nodes_per_strain, Metadata, by.x = "ATTRIBUTE_Strain", by.y = "ATTRIBUTE_Strain")
heatmap = heatmap %>% select(id, ATTRIBUTE_Strain, ATTRIBUTE_Bacteria)

## Create a list of unique IDs for each combination of strain and bacteria
df_grouped <- heatmap %>% 
  group_by(ATTRIBUTE_Strain, ATTRIBUTE_Bacteria) %>% 
  summarise(id_list = list(unique(id)), .groups = 'drop')

## Create a nested list of IDs by bacteria and strain
id_lists <- split(df_grouped$id_list, df_grouped$ATTRIBUTE_Bacteria)

## Obtain the names of the bacteria and the strains
bacteria_names <- names(id_lists)

## Create a list of intersection matrices for each bacteria
intersect_matrices <- list()

for (bacteria in bacteria_names) {
  id_per_strain <- id_lists[[bacteria]]
  strain_names <- df_grouped$ATTRIBUTE_Strain[df_grouped$ATTRIBUTE_Bacteria == bacteria]
  n_strains <- length(strain_names)
  
  # Intersection matrix for the specific bacteria
  intersect_matrix <- matrix(0, nrow = n_strains, ncol = n_strains, dimnames = list(strain_names, strain_names))
  
  for (i in 1:n_strains) {
    for (j in 1:n_strains) {
      intersect_matrix[i, j] <- length(intersect(id_per_strain[[i]], id_per_strain[[j]]))
    }
  }
  
  # Store the matrix in a list
  intersect_matrices[[bacteria]] <- intersect_matrix
}

## Convert matrices to long format
intersect_matrix_melt <- do.call(rbind, lapply(names(intersect_matrices), function(bacteria) {
  mat_melt <- as.data.frame(as.table(intersect_matrices[[bacteria]]))
  mat_melt$Bacteria <- bacteria
  return(mat_melt)
}))

## Plot the heatmap with facets by bacteria
Heatmap <- ggplot(intersect_matrix_melt, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("white", "lightblue", "blue", "darkblue")) +
  labs(x = "", y = "", fill = "Shared nodes number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        axis.text.y = element_text(size = 4)) 

Heatmap
## Save the Distribution Plot as a PNG file
ggsave("Data/IQNPA_Plot/4_Heatmap.png", plot = Heatmap, width = 15, height = 8)


# Network Analysis
library(igraph)
library(ggraph)
library(proxy)
library(visNetwork)

## Convert to matrix
data_matrix <- as.matrix(IQ_NPA[, c("Hmqs", "Hcos", "Ans")])

## Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

## Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ_NPA$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ_NPA$ATTRIBUTE_Strain

## Define threshold
threshold <- 0.995

## Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

## Convert to numeric matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

## Create the graph from the adjacency matrix
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

## Add ATTRIBUTE_Bacteria as a vertex attribute
V(g)$bacteria <- IQ_NPA$ATTRIBUTE_Bacteria

## Export the graph as a GraphML file
write_graph(g, file = "Data/IQNPA_Plot/Similarity.graphml", format = "graphml")

## Create a color vector based on ATTRIBUTE_Bacteria
IQ_NPA$ATTRIBUTE_Bacteria <- as.factor(IQ$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ$ATTRIBUTE_Bacteria)))[as.numeric(IQ$ATTRIBUTE_Bacteria)]

## Visualize the network
plot(g, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 5, 
     vertex.color = color_vector,  # Use the color vector
     main = "")


# Summary of top 10 lowest IQ
IQ_Summary_top <- IQ_NPA %>%
  arrange(IQ_NPA) %>%
  slice_head(n = 10)

## IQ top plot
IQ_top_plot = ggplot(IQ_Summary_top, aes(x = reorder(ATTRIBUTE_Strain, -IQ), y = IQ, fill = ATTRIBUTE_Bacteria)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "Black", width = 0.7, show.legend = TRUE) +  
  theme_classic() + 
  labs(x = "",
       y = "IQ",
       fill = "Bacteria") +  # Legend label
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 14),  
    axis.title = element_text(face = "bold"),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),  
    legend.title = element_text(size = 14)  
  )+ scale_fill_viridis_d()

IQ_top_plot
## Save the plot
ggsave("Data/IQNPA_Plot/5_IQ_top_Plot.png", plot = IQ_top_plot, width = 15, height = 8)

## Distribution top plot data
summary_top_data <- IQ_Summary_top %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(
    mean_IQ = mean(IQ, na.rm = TRUE),
    sd_IQ = sd(IQ, na.rm = TRUE),
    .groups = 'drop'
  )

## Distribution top plot
Distribution_top_plot = ggplot() +
  geom_point(data = IQ_Summary_top, aes(x = ATTRIBUTE_Bacteria, y = IQ), 
             shape = 1, color = "blue", size = 1.5, alpha = 0.6) +  # Points for individual IQ values
  geom_errorbar(data = summary_top_data, 
                aes(x = ATTRIBUTE_Bacteria, ymin = mean_IQ - sd_IQ, ymax = mean_IQ + sd_IQ), 
                width = 0.1, color = "Black") +  # Error bars for mean IQ
  geom_point(data = summary_top_data, aes(x = ATTRIBUTE_Bacteria, y = mean_IQ), 
             shape = 1, color = "blue", size = 1) +  # Points for mean IQ
  labs(
    x = "Bacteria", 
    y = "IQ") +
  theme_classic() 

Distribution_top_plot
# Save the distribution plot
ggsave("Data/IQNPA_Plot/6_Distribution_top_plot.png", plot = Distribution_top_plot, width = 15, height = 8)

# Heatmap top

heatmap_top = merge(heatmap, IQ_Summary_top, by.x = "ATTRIBUTE_strain", by.y = "ATTRIBUTE_Strain")
heatmap_top = heatmap_top %>% select(id, ATTRIBUTE_strain, ATTRIBUTE_Bacteria.x)
colnames(heatmap_top) = c("id", "ATTRIBUTE_strain", "ATTRIBUTE_Bacteria")

## Create a list of unique IDs for each strain and bacteria combination
df_grouped_top <- heatmap_top %>% 
  group_by(ATTRIBUTE_strain, ATTRIBUTE_Bacteria) %>% 
  summarise(id_list = list(unique(id)))

## Create a nested list of IDs by bacteria and strain
id_lists_top <- split(df_grouped_top$id_list, df_grouped_top$ATTRIBUTE_Bacteria)

## Get the names of the bacteria and strains
bacteria_names_top <- names(id_lists_top)

## Create a list of intersection matrices for each bacteria
intersect_matrices_top <- list()

for (bacteria in bacteria_names_top) {
  id_per_strain_top <- id_lists_top[[bacteria]]
  strain_names_top <- df_grouped_top$ATTRIBUTE_strain[df_grouped_top$ATTRIBUTE_Bacteria == bacteria]
  n_strains_top <- length(strain_names_top)
  
  # Intersection matrix for the specific bacteria
  intersect_matrix <- matrix(0, nrow = n_strains_top, ncol = n_strains_top, dimnames = list(strain_names_top, strain_names_top))
  for (i in 1:n_strains_top) {
    for (j in 1:n_strains_top) {
      intersect_matrix[i, j] <- length(intersect(id_per_strain_top[[i]], id_per_strain_top[[j]]))
    }
  }
  
  # Store the matrix in a list
  intersect_matrices_top[[bacteria]] <- intersect_matrix
}

## Convert the matrices to long format
intersect_matrix_melt_top <- do.call(rbind, lapply(names(intersect_matrices_top), function(bacteria) {
  mat_melt <- as.data.frame(as.table(intersect_matrices_top[[bacteria]]))  # Updated to intersect_matrices_top
  mat_melt$Bacteria <- bacteria
  return(mat_melt)
}))

## Plot the heatmap with facets by bacteria
Heatmap_top <- ggplot(intersect_matrix_melt_top, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colours = c("white", "lightblue", "blue", "darkblue")) +
  labs(x = "", y = "", fill = "Shared nodes number") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) 
Heatmap_top

# Save the distribution plot
ggsave("Data/IQNPA_Plot/7_Heatmap_top_plot.png", plot = Heatmap_top, width = 15, height = 8)


# Red analysis top 10

# Convert to matrix
data_matrix <- as.matrix(IQ_Summary_top[, c("Hmqs", "Hcos", "Ans")])

# Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

# Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ_Summary_top$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ_Summary_top$ATTRIBUTE_Strain

# Define threshold
threshold <- 0.995

# Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

# Convert to numerical matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

# Create the graph from the adjacency matrix
g <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Ensure ATTRIBUTE_Bacteria is character or factor
IQ_Summary_top$ATTRIBUTE_Bacteria <- as.character(IQ_Summary_top$ATTRIBUTE_Bacteria)

# Add ATTRIBUTE_Bacteria as a vertex attribute
V(g)$bacteria <- IQ_Summary_top$ATTRIBUTE_Bacteria

# Export the graph as a GraphML file
write_graph(g, file = "Data/IQNPA_Plot/Top_graphlm.graphml", format = "graphml")

# Create a color vector based on ATTRIBUTE_Bacteria
IQ_Summary_top$ATTRIBUTE_Bacteria <- as.factor(IQ_Summary_top$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ_Summary_top$ATTRIBUTE_Bacteria)))[as.numeric(IQ$ATTRIBUTE_Bacteria)]

# Visualize the network
plot(g, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 30, 
     vertex.color = color_vector,  # Use the color vector
     main = "")

# Proportions plot

# Count how many ATTRIBUTE_strain there are for each ATTRIBUTE_Bacteria in the original dataframe
strain_count_per_bacteria <- Metadata %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(n_strains = n_distinct(ATTRIBUTE_strain))

# Count how many ATTRIBUTE_strain there are for each ATTRIBUTE_Bacteria in the top dataframe
strain_count_per_bacteria_top <- IQ_Summary_top %>%
  group_by(ATTRIBUTE_Bacteria) %>%
  summarise(n_strains_top = n_distinct(ATTRIBUTE_Strain))

# Join both dataframes to calculate the proportion
proportion_counts <- strain_count_per_bacteria %>%
  left_join(strain_count_per_bacteria_top, by = "ATTRIBUTE_Bacteria") %>%
  mutate(n_strains_top = replace_na(n_strains_top, 0),  # Replace NA with 0
         proportion = n_strains_top / n_strains * 100)
proportion_counts$n_strain_top_proportion<- proportion_counts$n_strains_top/10*100

# Plot the proportions
PPlot<- ggplot(proportion_counts, aes(x = reorder(ATTRIBUTE_Bacteria, -proportion), y = proportion, fill = proportion)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "lightblue", high = "cadetblue") + 
  labs(
    x = "",
    y = "Strain (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8))
PPlot
ggsave("Data/IQNPA_Plot/8_Proportion_plot.png", plot = PPlot, width = 15, height = 8)


Bb_plot <- ggplot(proportion_counts, aes(x = proportion, y = n_strains_top, size = n_strains, color = ATTRIBUTE_Bacteria, label = ATTRIBUTE_Bacteria)) +
  geom_point(alpha = 0.7) + 
  geom_text(vjust = -2, hjust = 0.5, check_overlap = TRUE, size = 3) +  
  scale_size_continuous(range = c(3, 10)) +  
  scale_color_brewer(palette = "Set2") +  
  labs(x = "Relative proportion (%) [Number of strain in top / total strain]", 
       y = "Proportion in top (%) [Number of strain in top / top (n)]", 
       size = "Total strains") +  
  theme_classic() +  
  theme(legend.position = "right",  
        axis.title = element_text(size = 12),  
        axis.text = element_text(size = 10),  
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  guides(color = "none")  
Bb_plot
ggsave("Data/IQNPA_Plot/9_Proportion_Bubble_plot.png", plot = Bb_plot, width = 15, height = 8)


# Bini integration

# BiNi result

BiNI_results <- read_excel("Data/BiNI_results.xlsx")


IQ_BiNI <- merge(IQ_NPA, BiNI_results, by.x = "ATTRIBUTE_Strain", by.y = "Strain")
IQ_BiNI_long <- IQ_BiNI %>% select(ATTRIBUTE_Bacteria, ATTRIBUTE_Strain, IQ, BiNI)

# IQ-BiNi plot
IQ_BiNI_long <- IQ_BiNI_long %>%
  pivot_longer(
    cols = -c(ATTRIBUTE_Strain, ATTRIBUTE_Bacteria),
    names_to = "Parameter",    
    values_to = "Value"        
  )

# Create a custom cadet blue color gradient
palette <- colorRampPalette(c("#5F9EA0", "#4682B4", "#2F4F4F"))(length(unique(IQ_BiNI_long$ATTRIBUTE_Bacteria)))

# Create the plot with dynamically generated aquamarine colors
## Create the plot
IQ_BiNI_Plot <- ggplot(IQ_BiNI_long, aes(x = ATTRIBUTE_Strain, y = Value, color = ATTRIBUTE_Bacteria, shape = Parameter)) +
  geom_point(size = 3, alpha = 0.8) +  # Add transparency for better overlap handling
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_classic(base_size = 14) +  # Base font size increased for better readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),  
    axis.title.y = element_text(size = 12, margin = margin(r = 10))
  ) +
  ylim(0, 1) +
  scale_shape_manual(values = c(16, 17)) + 
  scale_color_manual(values = palette) +  
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
  scale_color_viridis_d() 

IQ_BiNI_Plot

ggsave("Data/IQNPA_Plot/10_IQ_BiNI_plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)


## Create the scatter plot of BiNI vs IQ for each strain
IQ_BiNI_Plot <- ggplot(IQ_BiNI, aes(x = BiNI, y = IQ, color = ATTRIBUTE_Bacteria, shape = ATTRIBUTE_Strain)) +
  geom_point(size = 3) +
  labs(title = "",
       x = "BiNI",
       y = "IQ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylim(0, 1) +
  xlim(0, 1) +  
  scale_shape_manual(values = 1:length(unique(IQ_BiNI$ATTRIBUTE_Strain))) +  
  scale_fill_viridis_d()  

IQ_BiNI_Plot

ggsave("Data/IQNPA_Plot/11_IQ_BiNI_s_Plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)

# scatter plot IQ_BiNI
# Define a color palette

palette <- brewer.pal(n = length(unique(IQ_BiNI$ATTRIBUTE_Bacteria)), name = "Set2")
correlation <- cor(IQ_BiNI$IQ, IQ_BiNI$BiNI, use = "complete.obs")
correlation

# Crear el gráfico
IQ_BiNI_SPlot <- ggplot(IQ_BiNI, aes(x = IQ, y = BiNI, color = ATTRIBUTE_Bacteria)) +
  geom_point(size = 4, alpha = 0.7) +  # Adjusting point transparency and size
  geom_smooth(method = "lm", linetype = "dashed", color = "black", fill = "lightgray") +  # Confidence interval
  labs(title = "",
       x = "IQ",
       y = "BiNI") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),  # Increase axis text size
    axis.title = element_text(size = 14),  # Increase axis title size
    plot.title = element_text(hjust = 0.5, size = 16)  # Center and increase plot title size
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("R: ", round(correlation, 2)), 
           hjust = 3, vjust = 10, size = 5, color = "black") +
  scale_color_viridis_d() 
IQ_BiNI_SPlot


IQ_BiNI_SPlot

ggsave("Data/IQNPA_Plot/11_IQ_BiNI_SPlot.png", plot = IQ_BiNI_SPlot, width = 15, height = 8)


# Convert to matrix
IQ_BiNI_f <- IQ_BiNI %>% select(-BiNI)
data_matrix <- as.matrix(IQ_BiNI_f[, c("Hmqs", "Hcos", "Ans")])

# Calculate the similarity matrix using cosine distance
similarity_matrix <- 1 - proxy::dist(data_matrix, method = "cosine")

# Convert to matrix and assign row and column names
similarity_matrix <- as.matrix(similarity_matrix)
rownames(similarity_matrix) <- IQ_BiNI_f$ATTRIBUTE_Strain
colnames(similarity_matrix) <- IQ_BiNI_f$ATTRIBUTE_Strain

# Define threshold
threshold <- 0.9

# Create the adjacency matrix based on the threshold
adjacency_matrix <- similarity_matrix > threshold  # TRUE/FALSE

# Convert to numerical matrix (0 and 1)
adjacency_matrix <- as.numeric(adjacency_matrix)
adjacency_matrix <- matrix(adjacency_matrix, nrow = nrow(similarity_matrix), ncol = ncol(similarity_matrix))
rownames(adjacency_matrix) <- rownames(similarity_matrix)
colnames(adjacency_matrix) <- colnames(similarity_matrix)

# Create the graph from the adjacency matrix
IQ_BiNI_G <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Ensure ATTRIBUTE_Bacteria is character or factor
IQ_BiNI_f$ATTRIBUTE_Bacteria <- as.character(IQ_BiNI_f$ATTRIBUTE_Bacteria)

# Add ATTRIBUTE_Bacteria as a vertex attribute
V(IQ_BiNI_G)$bacteria <- IQ_BiNI_f$ATTRIBUTE_Bacteria

# Export the graph as a GraphML file
write_graph(IQ_BiNI_G, file = "Data/IQNPA_Plot/IQ_BiNi_graphml.graphml", format = "graphml")

# Create a color vector based on ATTRIBUTE_Bacteria
IQ_BiNI_f$ATTRIBUTE_Bacteria <- as.factor(IQ_BiNI_f$ATTRIBUTE_Bacteria)
color_vector <- rainbow(length(levels(IQ_BiNI_f$ATTRIBUTE_Bacteria)))[as.numeric(IQ_BiNI_f$ATTRIBUTE_Bacteria)]

## Visualize the network
plot(IQ_BiNI_G, 
     vertex.label = V(g)$name, 
     edge.width = E(g)$weight * 5, 
     edge.color = "gray", 
     vertex.size = 30, 
     vertex.color = color_vector,  # Use the color vector
     main = "")

# IQ vs IQNP

# Combine the two data frames by matching on ATTRIBUTE_Strain and ATTRIBUTE_Bacteria
combined_data <- data.frame(
  ATTRIBUTE_Strain = IQ$ATTRIBUTE_Strain,
  ATTRIBUTE_Bacteria = IQ$ATTRIBUTE_Bacteria,
  IQ = IQ$IQ,
  IQ_NP = IQ_NPA$IQ
)
colnames(IQ_combined)=c("ATTRIBUTE_Bacteria", "IQ", "IQNP")
# Plotting IQ vs IQ_NP, differentiating by ATTRIBUTE_Bacteria

correlation <- cor(combined_data$IQ, combined_data$IQNP, use = "complete.obs")
correlation
IQ_vs_IQ_NP_Plot <- ggplot(combined_data, aes(x = IQ, y = IQNP, color = ATTRIBUTE_Bacteria)) +
  geom_point(size = 4, alpha = 0.7) +  # Adjusting point transparency and size
  geom_smooth(method = "lm", linetype = "dashed", color = "black", fill = "lightgray") +  # Confidence interval
  labs(title = "",
       x = "IQ",
       y = "IQNP") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),  
    axis.title = element_text(size = 14),  
    plot.title = element_text(hjust = 0.5, size = 16)  
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("R: ", round(correlation, 2)), 
           hjust = 3, vjust = 10, size = 5, color = "black") +
  scale_color_viridis_d() 
IQ_vs_IQ_NP_Plot
ggsave("Data/IQNPA_Plot/12_IQ_IQNP_SPlot.png", plot = IQ_vs_IQ_NP_Plot, width = 15, height = 8)


# parameters IQ plot
IQ_long <- IQ %>%
  pivot_longer(cols = c("Ans", "Hcos", "Hmqs"), 
               names_to = "Parameter", 
               values_to = "Value") %>%
  mutate(Value = Value / 3)  
IQ_long <- IQ_long %>%
  group_by(ATTRIBUTE_Strain) %>%
  mutate(Total = sum(Value)) %>%
  ungroup() %>%
  arrange(desc(Total))


IQ_plot <- ggplot(IQ_long, aes(x = reorder(ATTRIBUTE_Strain, -Total), y = Value, fill = Parameter)) +
  geom_bar(stat = "identity") +  # Stacked bars
  labs(x = "Strain", 
       y = "Parameter value", 
       title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),  
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 18) 
  ) +
  scale_fill_brewer(palette = "Set2")  


print(IQ_plot)
ggsave("Data/IQ_Plot/12_IQ_input_plot.png", plot = IQ_plot, width = 15, height = 8)


# parameters IQNP plot
IQ_NPA_long <- IQ_NPA %>%
  pivot_longer(cols = c("Ans", "Hcos", "Hmqs"), 
               names_to = "Parameter", 
               values_to = "Value")%>%
  mutate(Value = Value / 3)

IQ_NPA_long  <- IQ_NPA_long  %>%
  group_by(ATTRIBUTE_Strain) %>%
  mutate(Total = sum(Value)) %>%
  ungroup() %>%
  arrange(desc(Total))  

IQ_NPA_plot  <- ggplot(IQ_NPA_long , aes(x = reorder(ATTRIBUTE_Strain, -Total), y = Value, fill = Parameter)) +
  geom_bar(stat = "identity") +  # Stacked bars
  labs(x = "Strain", 
       y = "Input value", 
       title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 16),  # Adjust x-axis labels
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.text = element_text(size = 20),  
    legend.title = element_text(size = 18) 
  ) +
  scale_fill_brewer(palette = "Set2")  # Use a color palette for different parameters


print(IQ_NPA_plot)
ggsave("Data/IQNPA_Plot/12_IQNP_input_plot.png", plot = IQ_NPA_plot, width = 15, height = 8)


# IQ and IQ_NPS distribution

# Combine the two data frames into one for easier plotting
IQ_combined <- IQ %>%
  mutate(Source = "IQ") %>%  # Add a column to indicate the source
  select(ATTRIBUTE_Bacteria, IQ, Source) %>%
  bind_rows(IQ_NPA %>% mutate(Source = "IQNP") %>% select(ATTRIBUTE_Bacteria, IQ, Source))


# Calculate summary statistics for both IQ and IQ_NPA
summary_data <- IQ_combined %>%
  group_by(Source, ATTRIBUTE_Bacteria) %>%
  summarise(
    mean_IQ = mean(IQ, na.rm = TRUE),  # Calculate mean IQ
    sd_IQ = sd(IQ, na.rm = TRUE),      # Calculate standard deviation of IQ
    .groups = 'drop'                   # Drop grouping structure
  )

colnames(IQ_combined)<-c("ATTRIBUTE_Bacteria", IQ, Index)
# Distribution Plot

Distribution_plot = ggplot(IQ_combined, aes(x = ATTRIBUTE_Bacteria, y = IQ, fill = Source)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +  
  labs(
    x = "", 
    y = "",
    fill="Index") + 
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"), 
    axis.title.y = element_text(size = 16, face = "bold"),  
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 14), 
    axis.line = element_line(size = 1), 
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 12),  
    strip.text = element_text(size = 16, face = "bold"),  
    strip.background = element_rect(fill = "WHITE", color = "black"), 
    plot.title = element_text(hjust = 0.5, size = 18)  
  ) +
  scale_fill_manual(values = c("IQ" = "#008B8B", "IQNP" = "#A52A2A")) 

  # Print the plot
print(Distribution_plot)

ggsave("Data/IQNPA_Plot/14_IQ_IQNP_BOXPLOT_plot.png", plot = Distribution_plot, width = 15, height = 8)



