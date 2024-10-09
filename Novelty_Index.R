
#Packages instalation

if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse", "readxl", "rvest","dplyr","tidyr","igraph","visNetwork", "readr","ggplot2", "ggraph", "graphTweets", "xmls", "writexl")

Sys.getenv("JAVA_HOME")
Sys.setenv(JAVA_HOME='C:/Program Files/Java/jdk-17')
library(rJava)
library(rcdk)
library(RColorBrewer)
library(writexl)

# Load network from a GraphML file. Change name

graphml <-  read_graph("Data/network_Crusemann.graphml", format = "graphml")

# Get node and edge data
if (is_igraph(graphml)) {
  df_edges <- get.data.frame(graphml, what = "edges")
  df_nodes <- get.data.frame(graphml, what = "vertices")
}

#Extract index inputs
score_MSQ <- df_nodes %>%
  select(id, MQScore, UniqueFileSources) %>%
  filter(!is.na(MQScore) & MQScore != "") # Extract MSQ and remove NAs
score_COS_node1 <- df_edges %>%
  select(node1, cosine_score, node2) %>%
  semi_join(score_MSQ, by = c("node2" = "id")) #Cosine values inputs from node to node
score_COS_node2 <- df_edges %>%
  select(node2, cosine_score, node1) %>%
  semi_join(score_MSQ, by = c("node1" = "id")) #Complementary cosine values inputs

# Jump to line 500 for an index based on NPAtlas annotation nodes

# Obtain Cos value
combined_interactions <- bind_rows(
  score_COS_node1 %>% mutate(direction = "Score_COS_node1"),
  score_COS_node2 %>% mutate(direction = "Score_COS_node2")
) #Add a new column called direction to indicate the direction of the interaction by referencing the node of interest.

combined_interactions <- combined_interactions %>% mutate(
  combined_nodes = ifelse(direction == "Score_COS_node1",
                          paste(node1, node2, sep = "/"),
                          paste(node2, node1, sep = "/"))
) #Verification of the value of the direction column. If “Score_COS_node1”, concatenate node1 and node2, 
#concatenates node1 and node2. If “Score_COS_node2”, concatenates node2 and node1.

combined_interactions <- combined_interactions %>%
  separate(combined_nodes, into = c("nodeA", "nodeB"), sep = "/")# Now node A takes the position of the node of interest.


# Eliminate rows where nodeA and nodeB are the same, i.e. singlets
combined_interactions <- combined_interactions[combined_interactions$nodeA != combined_interactions$nodeB, ]

score_COS <- combined_interactions %>%
  select(nodeA, cosine_score) # The Cosine value is filtered according to the node of interest (attached to an annotated node).

# Merge dataframes by id
BASE <- merge(score_MSQ, score_COS, by.x = "id", by.y = "nodeA", all = TRUE) # Note that the MQS value is repeated when there is more than one COS, which is exactly what was sought


########################## Retrieving the third parameter, number of unannotated nodes
score_MSQ <- df_nodes %>%
  select(id, MQScore, UniqueFileSources) %>%
  filter(!is.na(MQScore) & MQScore != "") %>%
  select(-MQScore)

df_nodes_annotated <- score_MSQ %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>% # Change Sample identifier
  group_by(UniqueFileSources)

# Count id annotated per sample
count_nodes_annotated_sample <- table(df_nodes_annotated$UniqueFileSources)
df_count<- as.data.frame(count_nodes_annotated_sample)

# total nodes per sample
nodes_sample <- df_nodes %>%
  select(id, UniqueFileSources) %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>% # Change Sample identifier
  group_by(UniqueFileSources)

# Count total nodes per sample
count_nodes_sample <- table(nodes_sample$UniqueFileSources)
df_count_total_nodes<- as.data.frame(count_nodes_sample)

# Calculate nodes annotated per sample/ total nodes per sample
nodes_annotated_total= merge(df_count, df_count_total_nodes, by= "Var1")
nodes_annotated_total$Percentage <- (nodes_annotated_total$Freq.x / nodes_annotated_total$Freq.y) 
df_Unn <- nodes_annotated_total %>% select(-Freq.x,-Freq.y)

######Shannon index for MQS y Cos
# Verificar si el directorio existe, si no, crearlo
if (!dir.exists("Data/CSV_IQ")) {
  dir.create("Data/CSV_IQ")
}

# Separate the samples and extend the dataframe
df_separated <- BASE %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>%
  group_by(UniqueFileSources)

# Define strata for MQScore
MQS_breaks <- seq(0.01, 1, by = 0.05)
strata_MQS <- setNames(
  lapply(seq_along(MQS_breaks)[-length(MQS_breaks)], function(i) c(MQS_breaks[i], MQS_breaks[i + 1])),
  paste0(sprintf("%.2f", MQS_breaks[-length(MQS_breaks)]), "-", sprintf("%.2f", MQS_breaks[-1])))


# Define strata for cosine score
cos_breaks <- seq(0.70, 1.00, by = 0.05)
strata_COS <- setNames(
  lapply(seq_along(cos_breaks)[-length(cos_breaks)], function(i) c(cos_breaks[i], cos_breaks[i + 1])),
  paste0(sprintf("%.2f", cos_breaks[-length(cos_breaks)]), "-", sprintf("%.2f", cos_breaks[-1])))

# Function to stratify and count values, ignoring NAs
stratify_and_count <- function(df, column, strata) {
  counts <- sapply(strata, function(range) {
    sum(df[[column]] >= range[1] & df[[column]] < range[2], na.rm = TRUE)
  })
  return(counts)
}

# Function to stratify and count values, ignoring NAs
stratify_and_count_unique <- function(df, column, strata) {
  unique_df <- df %>% 
    group_by(id) %>% 
    slice(1) %>%  # Mantener la primera ocurrencia por id
    ungroup()
  
  # count unique values for strata
  counts <- sapply(strata, function(range) {
    sum(unique_df[[column]] >= range[1] & unique_df[[column]] < range[2], na.rm = TRUE)
  })
  return(counts)
}

# Function to calculate Shannon Index
calculate_shannon_index <- function(counts) {
  proportions <- counts / sum(counts)
  proportions <- proportions[proportions > 0] # Eliminate proportions 0 to avoid log(0)
  H <- -sum(proportions * log(proportions))
  return(H)
}

# Create lists to save individual df´s and strata counts
individual_dfs <- list()
strata_counts <- list()

# Create a df´s to save Shannon Index for each sample
IQ_Summary <- data.frame(
  Sample = character(),
  MQScore_Shannon_Index = numeric(),
  Cosine_Score_Shannon_Index = numeric(),
  stringsAsFactors = FALSE
)

# For each sample, create an individual df and strata counts
for (sample_name in unique(df_separated$UniqueFileSources)) {
  individual_df <- df_separated %>%
    filter(UniqueFileSources == sample_name) %>%
    select(id, MQScore, cosine_score)
  
  # Counting unique values in each stratum for MQScore
  MQScore_counts <- stratify_and_count_unique(individual_df, "MQScore", strata_MQS)
  
  # Count values in each stratum for cosine_score
  cosine_score_counts <- stratify_and_count(individual_df, "cosine_score", strata_COS)
  
  # Calculating the Shannon index for MQScore
  MQScore_shannon <- calculate_shannon_index(MQScore_counts)
  
  # Calculate Shannon index for cosine_score
  cosine_score_shannon <- calculate_shannon_index(cosine_score_counts)
  
  # Ensure that the strata names are consistent.
  counts_df <- data.frame(
    Stratum = names(strata_MQS),
    MQScore_Count = MQScore_counts,
    Cosine_Score_Count = cosine_score_counts[1:length(MQScore_counts)]  
  )
  
  # Añadir el índice de Shannon como la última fila
  counts_df <- rbind(counts_df, data.frame(
    Stratum = "Shannon_Index",
    MQScore_Count = MQScore_shannon,
    Cosine_Score_Count = cosine_score_shannon
  ))
  
  # Store count dataframe in the list
  strata_counts[[sample_name]] <- counts_df
  
  # Assign the individual dataframe to the list
  individual_dfs[[sample_name]] <- individual_df
  
  # Adding the Shannon index to the summary dataframe
  IQ_Summary <- rbind(IQ_Summary, data.frame(
    Sample = sample_name,
    MQScore_Shannon_Index = MQScore_shannon,
    Cosine_Score_Shannon_Index = cosine_score_shannon
  ))
}


# Merge the IQ_Summary with df_conunt based on the Sample and Var1 columns
IQ_Summary <- merge(IQ_Summary, df_Unn, by.x = "Sample", by.y = "Var1", all.x = TRUE)

# Rename columns and reorder them
colnames(IQ_Summary) <- c("Sample", "Hmsq", "Hcos", "Ans")
IQ_Summary <- IQ_Summary[, c("Sample", "Hmsq", "Hcos", "Ans")]

# Delete rows with NA in columns Hmsq and Hcos
IQ_Summary <- IQ_Summary[complete.cases(IQ_Summary$Hmsq, IQ_Summary$Hcos), ]


# min-max normalization
IQ_Summary$Hmsq <- (IQ_Summary$Hmsq - min(IQ_Summary$Hmsq)) / (max(IQ_Summary$Hmsq) - min(IQ_Summary$Hmsq))
IQ_Summary$Hcos <- (IQ_Summary$Hcos - min(IQ_Summary$Hcos)) / (max(IQ_Summary$Hcos) - min(IQ_Summary$Hcos))


IQ_Summary$IQ <- ifelse(is.na(IQ_Summary$Hmsq), 0, IQ_Summary$Hmsq) +
  ifelse(is.na(IQ_Summary$Hcos), 0, IQ_Summary$Hcos) +
  ifelse(is.na(IQ_Summary$Ans), 0, IQ_Summary$Ans)

IQ_Summary$IQ= (3-IQ_Summary$IQ)/3

# Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder (Se guardan con un nombre distinto)
for(sample_name in names(individual_dfs)) {
  write.csv(individual_dfs[[sample_name]], file = paste0("Data/CSV_IQ/", sample_name, "_data.csv"), row.names = FALSE)
  write.csv(strata_counts[[sample_name]], file = paste0("Data/CSV_IQ/", sample_name, "_strata_counts.csv"), row.names = FALSE)
  
  # Assign each dataframe to a variable in the R environment
  assign(sample_name, individual_dfs[[sample_name]])
}

# Save the updated Shannon Index summary dataframe to a CSV file in the CSV_IQ folder
write.csv(IQ_Summary, file = "Data/CSV_IQ/1_IQ_Summary.csv", row.names = FALSE)


# IQ plot for each sample

# Define the new number of samples per plot
samples_per_plot <- 10   # change it as required

# Recalculate the number of plots needed based on the new value
num_plots <- ceiling(nrow(IQ_Summary) / samples_per_plot)

if (!dir.exists("Data/IQ_Plot")) {
  dir.create("Data/IQ_Plot")
}

# Plots
for (i in 1:num_plots) {
  # Define the start and end index for the subset of data
  start_index <- (i - 1) * samples_per_plot + 1
  end_index <- min(i * samples_per_plot, nrow(IQ_Summary))
  # Subset the data for the current plot
  subset_data <- IQ_Summary[start_index:end_index, ]
  
  # Replace NA values with 0 in columns
  subset_data$Hmsq[is.na(subset_data$Hmsq)] <- 0
  subset_data$Hcos[is.na(subset_data$Hcos)] <- 0
  subset_data$Ans[is.na(subset_data$Ans)] <- 0
  
  # Create the plot
  p <- ggplot(subset_data, aes(x = Sample)) +
    geom_point(aes(y = Hmsq, color = "Hmsq"), size = 3) +
    geom_point(aes(y = Hcos, color = "Hcos"), size = 3, shape = 17) +
    geom_point(aes(y = Ans, color = "Ans"), size = 3, shape = 15) +
    labs(title = paste("IQ Index by Sample", start_index, "to", end_index),
         x = "Sample",
         y = "IQ Index",
         color = "Parameters") +
    theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)) +
    scale_color_manual(values = c("Hmsq" = "blue", "Hcos" = "red", "Ans" = "green")) +
    theme_classic() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  
  # Save the plot to a file in the IQ_Plot folder
  ggsave(paste0("Data/IQ_Plot/IQ_Plot_", i, ".png"), plot = p, width = 15, height = 8)
}


#Clustering plot

IQ_by_sample=IQ_Summary %>% select(-Hmsq, -Hcos, -Ans)

Metadata_Cruseman <- read_delim("Data/Metadata_Cruseman.csv", 
                                delim = ";", escape_double = FALSE, trim_ws = TRUE)

IQ_w_metadata<- merge(IQ_by_sample, Metadata_Cruseman, by.x = "Sample", by.y = "filename", all = TRUE)

IQ_w_metadata <- na.omit(IQ_w_metadata)

# Escala las variables de interés (IQ)
IQ_w_metadata_vector<- scale(IQ_w_metadata$IQ)

# Apply k-means with 3 cluster (modify as required)
set.seed(123)  # to reproducibility
kmeans_result <- kmeans(IQ_w_metadata_vector, centers = 10)

# Añadir la columna de clusters al dataframe
IQ_w_metadata$Cluster <- as.factor(kmeans_result$cluster)

# Plot for clusters. Modify ATTRIBUTES
IQ_Cluster_plot= ggplot(IQ_w_metadata, aes(x = Cluster, y = IQ, color = ATTRIBUTE_solvent, shape = ATTRIBUTE_Bacteria)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +  # Jitter en el eje x
  labs(title = "IQ by clusters with samples",
       x = "Clusters",
       y = "IQ",  
       color = "ATTRIBUTE_solvent",
       shape = "ATTRIBUTE_Bacteria") +  # Etiqueta para el shape
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set1")

# Save the plot to a file in the IQ_Plot folder
ggsave("Data/IQ_Plot/1_IQ_Cluster_plot.png", plot = IQ_Cluster_plot, width = 15, height = 8)


# Summary top 10 lower SI

IQ_Summary_top <- IQ_Summary %>%
  arrange(desc(IQ)) %>%
  slice_head(n = 10) 

# Summary plot inputs
p_summary <- ggplot(IQ_Summary_top, aes(x = Sample)) +
  geom_point(aes(y = Hmsq, color = "Hmsq"), size = 3) +
  geom_point(aes(y = Hcos, color = "Hcos"), size = 3, shape = 17) +
  geom_point(aes(y = Ans, color = "Ans"), size = 3, shape = 15) +
  labs(title = "Top 10 Samples with the Lowest IQ parameters",
       x = "Sample",
       y = "Value",
       color = "Parameters") +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("Hmsq" = "blue", "Hcos" = "red", "Ans" = "green")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave("Data/IQ_Plot/2_Parameters_summary.png", plot = p_summary, width = 15, height = 8)


# IQ Summary plot
p_summary <- ggplot(IQ_Summary_top, aes(x = Sample)) +
  geom_point(aes(y = IQ, color = "IQ"), size = 3) +
  labs(title = "Top 10 Samples with the Lowest IQ",
       x = "Sample",
       y = "IQ") +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave("Data/IQ_Plot/3_IQ_Plot_summary.png", plot = p_summary, width = 15, height = 8)


# BiNI integration

# Samples with BiNI
Samples_BiNI <- read_excel("Data/Samples_BiNI.xlsx")

# Filter the data based on the sample names vector
IQ_filtrated <- IQ_Summary %>%
  filter(Sample %in% Samples_BiNI$filename) %>% 
  select(-Hmsq, -Hcos, -Ans)
IQ_filtrated = merge(IQ_filtrated, Metadata_Cruseman, by.x = "Sample", by.y = "filename")

# BiNi result

BiNI_results <- read_excel("Data/BiNI_results.xlsx")


IQ_BiNI <- merge(IQ_filtrated, BiNI_results, by.x = "ATTRIBUTE_strain", by.y = "Strain")

# IQ-BiNi plot

IQ_BiNI_Plot <- ggplot(IQ_BiNI, aes(x = ATTRIBUTE_solvent, y = IQ, color = ATTRIBUTE_solvent, shape = ATTRIBUTE_Bacteria)) +
  geom_point(size = 3) +
  facet_wrap(~ ATTRIBUTE_strain, scales = "free_y") +
  geom_hline(data = IQ_BiNI, aes(yintercept = BiNI), linetype = "dashed", color = "red") +
  labs(title = "IQ by Strain and Solvent with BiNI value in line",
       x = "Solvent",
       y = "IQ") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylim(0, 1)

print(IQ_BiNI_Plot)

ggsave("Data/IQ_Plot/4_IQ_BiNI_Plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)
write_xlsx(IQ_BiNI, path = "Data/CSV_IQ/2_IQ_BINI.xlsx")

# IQ by strain

IQ_w_metadata <- IQ_w_metadata %>%
  select(-ATTRIBUTE_solvent, -Cluster)
IQ_by_strain <- IQ_w_metadata %>%
  group_by(ATTRIBUTE_strain, ATTRIBUTE_Bacteria) %>%
  summarise(total_IQ = sum(IQ, na.rm = TRUE)) 

# Save dataframe 
write.csv(IQ_by_strain, file = "Data/CSV_IQ/3_IQ_by_strain.csv", row.names = FALSE)

#plot IQ by strain

# IQ by strain. Modify ATTRIBUTES
IQ_by_strain_plot= ggplot(IQ_by_strain, aes(x = ATTRIBUTE_strain, y = total_IQ, color = ATTRIBUTE_Bacteria)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +  # Jitter en el eje x
  labs(title = "IQ by strain",
       x = "Strain",
       y = "IQ",  
       color = "ATTRIBUTE_Bacteria") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set1")

IQ_by_strain_plot

ggsave("Data/IQ_Plot/5_IQ_by_strain_Plot.png", plot = IQ_by_strain_plot, width = 15, height = 8)

# Clustering plot by strain

IQ_by_strain <- na.omit(IQ_by_strain)

# Scaling
IQ_by_strain_vector<- scale(IQ_by_strain$total_IQ)

# Apply k-means with 3 cluster (modify as required)
set.seed(123)  # to reproducibility
kmeans_result_strain <- kmeans(IQ_by_strain_vector, centers = 10)

# Añadir la columna de clusters al dataframe
IQ_by_strain$Cluster <- as.factor(kmeans_result_strain$cluster)

# Plot for clusters. Modify ATTRIBUTES
IQ_by_strain_Cluster_plot= ggplot(IQ_by_strain, aes(x = Cluster, y = total_IQ, color = ATTRIBUTE_Bacteria)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +  # Jitter en el eje x
  labs(title = "IQ by clusters with samples",
       x = "Clusters",
       y = "IQ",  
       color = "ATTRIBUTE_Bacteria") + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set1")

IQ_by_strain_Cluster_plot

# Save the plot to a file in the IQ_Plot folder
ggsave("Data/IQ_Plot/6_IQ_by_strain_Cluster_plot.png", plot = IQ_Cluster_plot, width = 15, height = 8)












# Summary top 10 lower total_IQ_by_strain

print(IQ_by_strain)

IQ_by_strain_Summary_top <- IQ_by_strain %>%
  arrange(desc(total_IQ)) %>%
  slice_head(n = 10) 


# IQ_by_strain Summary plot
IQ_by_strain_Summary_top_plot <- ggplot(IQ_by_strain_Summary_top, aes(x = ATTRIBUTE_strain)) +
  geom_point(aes(y = total_IQ, color = "IQ"), size = 3) +
  labs(title = "Top 10 Samples with the Lowest IQ by strain",
       x = "Sample",
       y = "IQ") +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.box = "horizontal")

IQ_by_strain_Summary_top_plot

ggsave("Data/IQ_Plot/3_IQ_Plot_summary.png", plot = p_summary, width = 15, height = 8)




1















#################################################
#INTEGRACION DE LOS INDICES POR CEPA
IQ_for_strain=IQ_Summary_NPAtlas %>% select(-Hmsq, -Hcos, -Nsna)

Metadata_Cruseman <- read_delim("Metadata_Cruseman.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)


IQ_NPAtlas_by_strain_w_metadata<- merge(IQ_for_strain, Metadata_Cruseman, by.x = "Sample", by.y = "filename", all = TRUE)
IQ_NPAtlas_by_strain_w_metadata<- aggregate(IQ ~ ATTRIBUTE_strain, data = IQ_NPAtlas_by_strain_w_metadata, FUN = sum)

Metadata_Cruseman=Metadata_Cruseman %>% select(-filename, -ATTRIBUTE_solvent)
IQ_NPAtlas_by_strain_w_metadata<- merge(IQ_NPAtlas_by_strain_w_metadata, Metadata_Cruseman, by.x = "ATTRIBUTE_strain", by.y =  "ATTRIBUTE_strain")
IQ_NPAtlas_by_strain_w_metadata_clean <- IQ_NPAtlas_by_strain_w_metadata %>% distinct()

# Clustering plot

IQ_NPAtlas_by_strain_w_metadata_clean <- na.omit(IQ_NPAtlas_by_strain_w_metadata_clean)

# Escala las variables de interés (IQ)
IQ_NPAtlas_by_strain_w_metadata_clean_scaled <- scale(IQ_NPAtlas_by_strain_w_metadata_clean$IQ)

# Aplicar k-means con 3 clusters (puedes ajustar el número de clusters)
set.seed(123)  # Para reproducibilidad
kmeans_result_NPAtlas_by_strain <- kmeans(IQ_NPAtlas_by_strain_w_metadata_clean_scaled, centers = 5)

# Añadir la columna de clusters al dataframe
IQ_NPAtlas_by_strain_w_metadata_clean$Cluster <- as.factor(kmeans_result_NPAtlas_by_strain$cluster)


library(RColorBrewer)

# Usar una paleta con colores diferenciados
ggplot(IQ_NPAtlas_by_strain_w_metadata_clean, aes(x = Cluster, y = IQ, color = ATTRIBUTE_Bacteria)) +
  geom_jitter(size = 3, width = 0.2, height = 0, alpha = 0.8) +  # Jitter en el eje x
  labs(title = "Visualización de IQ por Clusters",
       x = "Clusters",
       y = "IQ",  
       color = "Bacteria") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_brewer(palette = "Set1")  # Usar paleta Set1 de RColorBrewer















top_10_bajos <- IQ_for_strain %>% 
  arrange(IQ) %>%          
  slice_head(n = 10)

# Gráfico de puntos de IQ por ATTRIBUTE_strain
ggplot(top_10_bajos, aes(x = ATTRIBUTE_strain, y = IQ)) +
  geom_point(color = "blue", size = 3) +
  labs(title = "Top 10 lowest IQ by strain",
       x = "Strain (ATTRIBUTE_strain)",
       y = "IQ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje x






