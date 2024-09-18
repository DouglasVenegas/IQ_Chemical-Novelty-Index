if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse","dplyr","tidyr","igraph","visNetwork", "readr","ggplot2", "ggraph", "graphTweets", "xmls", "writexl")


# Load network from a GraphML file
graphml <-  read_graph("network_Crusemann.graphml", format = "graphml")

# Verify if the object is valid
if (is_igraph(graphml)) {
  # Get node and edge data
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


# Retrieving the third parameter, number of unannotated nodes
score_no_MSQ <- df_nodes %>%
  select(id, MQScore, UniqueFileSources) %>%
  filter(is.na(MQScore) | MQScore == "")%>%
  select(-MQScore)

df_nodes_unannotated <- score_no_MSQ %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>% # Change Sample identifier
  group_by(UniqueFileSources)

# Count id per sample
count <- table(df_nodes_unannotated$UniqueFileSources)
df_count<- as.data.frame(conteo)
total_freq <- sum(df_count$Freq)

# Calculate the percentage of each frequency
df_count$Percentage <- (df_count$Freq / total_freq)*100 # is maintained in percentage by the scale with the other parameters

df_count <- df_count %>% select(-Freq)

######Shannon index for MQS y Cos

if (!dir.exists("CSV_IQ")) {
  dir.create("CSV_IQ")
}


# Merge dataframes by id
BASE <- merge(score_MSQ, score_COS, by.x = "id", by.y = "nodeA", all = TRUE)

# Split samples and extend the dataframe
df_separated <- BASE %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>%
  group_by(UniqueFileSources)


# Define the strata

MQS_breaks <- seq(0.01, 1, by = 0.05)
strata_MQS <- setNames(
  lapply(seq_along(MQS_breaks)[-length(MQS_breaks)], function(i) c(MQS_breaks[i], MQS_breaks[i+1])),
  paste0(sprintf("%.2f", MQS_breaks[-length(MQS_breaks)]), "-", sprintf("%.2f", MQS_breaks[-1]))
)

# Define the strata for COS
cos_breaks <- seq(0.70, 1.00, by = 0.05)
strata_COS <- setNames(
  lapply(seq_along(cos_breaks)[-length(cos_breaks)], function(i) c(cos_breaks[i], cos_breaks[i+1])),
  paste0(sprintf("%.2f", cos_breaks[-length(cos_breaks)]), "-", sprintf("%.2f", cos_breaks[-1]))
)

# Function to stratify and count values, ignoring NAs
stratify_and_count <- function(df, column, strata) {
  counts <- sapply(strata, function(range) {
    sum(df[[column]] >= range[1] & df[[column]] < range[2], na.rm = TRUE)
  })
  return(counts)
}

# Function to stratify and count unique values, ignoring NAs
stratify_and_count_unique <- function(df, column, strata) {
  # Eliminate duplicate ids by keeping only one value per id
  unique_df <- df %>% 
    group_by(id) %>% 
    slice(1) %>%  # Keep the first occurrence for each id
    ungroup()
  # Count unique values in each strata
  counts <- sapply(strata, function(range) {
    sum(unique_df[[column]] >= range[1] & unique_df[[column]] < range[2], na.rm = TRUE)
  })
  return(counts)
}

# Function to calculate Shannon Index
calculate_shannon_index <- function(counts) {
  proportions <- counts / sum(counts)
  proportions <- proportions[proportions > 0] # Eliminate zero proportions to avoid NaN from log(0)
  H <- -sum(proportions * log(proportions))
  return(H)
}


# Create lists to store individual dataframes and strata counts
individual_dfs <- list()
strata_counts <- list()

# Create a dataframe to store the Shannon Index for each sample
IQ_Summary <- data.frame(
  Sample = character(),
  MQScore_Shannon_Index = character(),
  Cosine_Score_Shannon_Index = character(),
  stringsAsFactors = FALSE
)

# For each sample, create an individual dataframe and count the values in each strata
for(sample_name in unique(df_separated$UniqueFileSources)) {
  individual_df <- df_separated %>%
    filter(UniqueFileSources == sample_name) %>%
    select(id, MQScore, cosine_score)
  
  # Count values in each strata for MQScore
  MQScore_counts <- stratify_and_count_unique(individual_df, "MQScore", strata_MQS)
  
  # Count values in each strata for cosine_score
  cosine_score_counts <- stratify_and_count(individual_df, "cosine_score", strata_COS)
  
  # Calculate Shannon Index for MQScore
  MQScore_shannon <- calculate_shannon_index(MQScore_counts)
  print(MQScore_shannon)
  
  # Calculate Shannon Index for cosine_score
  cosine_score_shannon <- calculate_shannon_index(cosine_score_counts)
  
  # Combine the counts into a single dataframe
  counts_df <- data.frame(
    Stratum = names(strata),
    MQScore_Count = MQScore_counts,
    Cosine_Score_Count = cosine_score_counts
  )
  
  # Add the Shannon Index as the last row
  counts_df <- rbind(counts_df, data.frame(
    Stratum = "Shannon_Index",
    MQScore_Count = MQScore_shannon,
    Cosine_Score_Count = cosine_score_shannon
  ))
  
  # Store the counts dataframe in the list
  strata_counts[[sample_name]] <- counts_df
  
  # Assign the individual dataframe to the list
  individual_dfs[[sample_name]] <- individual_df
  
  # Add Shannon Index to the summary dataframe
  IQ_Summary <- rbind(IQ_Summary, data.frame(
    Sample = sample_name,
    MQScore_Shannon_Index = MQScore_shannon,
    Cosine_Score_Shannon_Index = cosine_score_shannon
  ))
}


# Merge the IQ_Summary with df_conteo based on the Sample and Var1 columns
IQ_Summary <- merge(IQ_Summary, df_conteo, by.x = "Sample", by.y = "Var1", all.x = TRUE)

# Rename columns and reorder them
colnames(IQ_Summary) <- c("Sample", "Hmsq", "Hcos", "Nsna")
IQ_Summary <- IQ_Summary[, c("Sample", "Hmsq", "Hcos", "Nsna")]

IQ_Summary$IQ <- ifelse(is.na(IQ_Summary$Hmsq), 0, IQ_Summary$Hmsq) +
  ifelse(is.na(IQ_Summary$Hcos), 0, IQ_Summary$Hcos) -
  ifelse(is.na(IQ_Summary$Nsna), 0, IQ_Summary$Nsna)

# Save each individual dataframe and the strata counts to CSV files in the CSV_IQ folder (Se guardan con un nombre distinto)
for(sample_name in names(individual_dfs)) {
  write.csv(individual_dfs[[sample_name]], file = paste0("CSV_IQ/", sample_name, "_data.csv"), row.names = FALSE)
  write.csv(strata_counts[[sample_name]], file = paste0("CSV_IQ/", sample_name, "_strata_counts.csv"), row.names = FALSE)
  
  # Assign each dataframe to a variable in the R environment
  assign(sample_name, individual_dfs[[sample_name]])
}

# Save the updated Shannon Index summary dataframe to a CSV file in the CSV_IQ folder
write.csv(IQ_Summary, file = "CSV_IQ/1_IQ_Summary.csv", row.names = FALSE)

# Define the new number of samples per plot
samples_per_plot <- 10   # Cambiar este valor según sea necesario

# Recalculate the number of plots needed based on the new value
num_plots <- ceiling(nrow(IQ_Summary) / samples_per_plot)


if (!dir.exists("IQ_Plot")) {
  dir.create("IQ_Plot")
}

# Plots
for (i in 1:num_plots) {
  # Subset the data for the current plot
  start_index <- (i - 1) * samples_per_plot + 1
  end_index <- min(i * samples_per_plot, nrow(IQ_Summary))
  subset_data <- IQ_Summary[start_index:end_index, ]
  
  # Replace NA values with 0 in Shannon index columns
  subset_data$Hmsq[is.na(subset_data$Hmsq)] <- 0
  subset_data$Hcos[is.na(subset_data$Hcos)] <- 0
  subset_data$Nsna[is.na(subset_data$Nsna)] <- 0
  
  # Create the plot
  p <- ggplot(subset_data, aes(x = Sample)) +
    geom_point(aes(y = Hmsq, color = "Hmsq"), size = 3) +
    geom_point(aes(y = Hcos, color = "Hcos"), size = 3, shape = 17) +
    geom_point(aes(y = Nsna, color = "Nsna"), size = 3, shape = 15) +
    labs(title = paste("IQ Index for Samples", start_index, "to", end_index),
         x = "Sample",
         y = "IQ Index",
         color = "Index Type") +
    theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)) +
    scale_color_manual(values = c("Hmsq" = "blue", "Hcos" = "red", "Nsna" = "green")) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  
  # Save the plot to a file in the IQ_Plot folder
  ggsave(paste0("IQ_Plot/IQ_Plot_", i, ".png"), plot = p, width = 15, height = 8)
}

# Summary top 10 lower SI
IQ_Summary_top <- IQ_Summary %>%
  mutate(IQ= Hmsq + Hcos - Nsna) %>%
  arrange(IQ) %>%
  head(10)

# Summary plot inputs
p_summary <- ggplot(IQ_Summary_top, aes(x = Sample)) +
  geom_point(aes(y = Hmsq, color = "Hmsq"), size = 3) +
  geom_point(aes(y = Hcos, color = "Hcos"), size = 3, shape = 17) +
  geom_point(aes(y = Nsna, color = "Nsna"), size = 3, shape = 15) +
  labs(title = "Top 10 Samples with the Lowest IQ parameters",
       x = "Sample",
       y = "IQ parameters",
       color = "Index Type") +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("Hmsq" = "blue", "Hcos" = "red", "Nsna" = "green")) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave("IQ_Plot/Parameters_summary.png", plot = p_summary, width = 15, height = 8)


# IQ Summary plot
p_summary <- ggplot(IQ_Summary_top, aes(x = Sample)) +
  geom_point(aes(y = IQ, color = "IQ"), size = 3) +
  labs(title = "Top 10 Samples with the Lowest IQ Index",
       x = "Sample",
       y = "IQ Index") +
  theme(axis.text.x = element_text(angle = -45, hjust = 1, vjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave("IQ_Plot/IQ_Plot_summary.png", plot = p_summary, width = 15, height = 8)






# Define the sample names vector and append the suffix
vector_nombres <- c(
  "T150-3-Me.mzXML", "T150-3-Bu.mzXML", "T150-3-EA.mzXML",
  "H643-3-Me.mzXML", "H643-3-Bu.mzXML", "H643-3-EA.mzXML",
  "S991-3-Me.mzXML", "S991-3-Bu.mzXML", "S991-3-EA.mzXML",
  "S863-3-Me.mzXML", "S863-3-Bu.mzXML", "S863-3-EA.mzXML",
  "T133-3-Me.mzXML", "T133-3-Bu.mzXML", "T133-3-EA.mzXML",
  "B440-3-Me.mzXML", "B440-3-Bu.mzXML", "B440-3-EA.mzXML",
  "T138-3-Me.mzXML", "T138-3-Bu.mzXML", "T138-3-EA.mzXML",
  "S960-3-Me.mzXML", "S960-3-Bu.mzXML", "S960-3-EA.mzXML",
  "T569-3-Me.mzXML", "T569-3-Bu.mzXML", "T569-3-EA.mzXML",
  "S744-3-Me.mzXML", "S744-3-Bu.mzXML", "S744-3-EA.mzXML",
  "S051-3-Me.mzXML", "S051-3-Bu.mzXML", "S051-3-EA.mzXML",
  "T250-3-Me.mzXML", "T250-3-Bu.mzXML", "T250-3-EA.mzXML",
  "S673-3-Me.mzXML", "S673-3-Bu.mzXML", "S673-3-EA.mzXML",
  "Y011-3-Me.mzXML", "Y011-3-Bu.mzXML", "Y011-3-EA.mzXML",
  "Y202-3-Me.mzXML", "Y202-3-Bu.mzXML", "Y202-3-EA.mzXML",
  "T799-3-Me.mzXML", "T799-3-Bu.mzXML", "T799-3-EA.mzXML",
  "R114-3-Me.mzXML", "R114-3-Bu.mzXML", "R114-3-EA.mzXML",
  "H877-3-Me.mzXML", "H877-3-Bu.mzXML", "H877-3-EA.mzXML",
  "T859-3-Me.mzXML", "T859-3-Bu.mzXML", "T859-3-EA.mzXML"
)

# Filter the data based on the sample names vector
IQ_filtrado <- IQ_Summary %>%
  filter(Sample %in% vector_nombres)

IQ_filtrado <- IQ_filtrado %>%
  mutate(Strain = sub("-.*", "", Sample))

IQ_filtrado <- IQ_filtrado %>%
  separate(Sample, into = c("Strain", "replica", "solvent_mzXML"), sep = "-", remove = FALSE) %>%
  separate(solvent_mzXML, into = c("solvent", "extension"), sep = "\\.", remove = FALSE) %>%
  select(-extension, -replica, -solvent_mzXML) 


# BiNi result
BiNi <- data.frame(
  Strain = c("T150", "H643", "S991", "S863", "T133", "B440", "T138", "S960", 
             "T569", "S744", "S051", "T250", "S673", "Y011", "Y202", "T799", 
             "R114", "H877", "T859"),
  BiNI = c(1.722, 1.653, 1.634, 1.701, 1.663, 1.773, 1.662, 1.854, 
           1.764, 1.677, 1.769, 1.863, 1.743, 1.579, 1.645, 1.590, 
           1.610, 1.694, 1.649)
)

df_combined <- merge(IQ_filtrado, BiNi, by = "Strain")

biNI_values <- df_combined %>%
  group_by(Strain) %>%
  summarise(BiNI_value = mean(BiNI, na.rm = TRUE))

# IQ-BiNi plot
IQ_BiNI_Plot<-ggplot(df_combined, aes(x = solvent, y = IQ, color = solvent)) +
  geom_point(size = 3) +
  facet_wrap(~ Strain, scales = "free_y") +
  geom_hline(data = biNI_values, aes(yintercept = BiNI_value), linetype = "dashed", color = "red") +
  labs(title = "IQ por Strain y Solvent con Línea de BiNI",
       x = "Solvent",
       y = "IQ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  ylim(1, 4) 

ggsave("IQ_Plot/IQ_BiNI_Plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)
write_xlsx(df_combined, path = "CSV_IQ/2_IQ_BINI.xlsx")


