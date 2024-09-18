library(dplyr)
library(tidyr)
library(igraph)
library(visNetwork)
library(readr)
library(ggplot2)
library(ggraph)
library(graphTweets)
library(xml2)
library(writexl)

#### Este script está adaptado a los datos de reanalisis del paper recomendado por parra
### la extracción de datos por cepa (en este caso) utiliza un comando específico para este .graphML, 
### se debe adaptar para otros muestras, si la muestra tiene un atributo único entonces se puede utilizar la metadata
### pero hay que pensar como adapatarlo de forma general


# Load network from a GraphML file
graphml <-  read_graph("network.graphml", format = "graphml")

# Verify if the object is valid
if (is_igraph(graphml)) {
  # Get node and edge data
  df_edges <- get.data.frame(graphml, what = "edges")
  df_nodes <- get.data.frame(graphml, what = "vertices")
}

# Extract index values

score_MSQ <- df_nodes %>%
  select(id, MQScore, UniqueFileSources) %>%
  filter(!is.na(MQScore) & MQScore != "") # Extract MSQ and remove NAs
score_COS_node1 <- df_edges %>%
  select(node1, cosine_score, node2) %>%
  semi_join(score_MSQ, by = c("node2" = "id")) #Valores coseno de un nodo a otro
score_COS_node2 <- df_edges %>%
  select(node2, cosine_score, node1) %>%
  semi_join(score_MSQ, by = c("node1" = "id")) #Valores coseno complementarios

# Obtener valor Cos
combined_interactions <- bind_rows(
  score_COS_node1 %>% mutate(direction = "Score_COS_node1"),
  score_COS_node2 %>% mutate(direction = "Score_COS_node2")
) #Añaden una nueva columna llamada direction para indicar la dirección de la interacción referenciando el nodo de interés

combined_interactions <- combined_interactions %>% mutate(
  combined_nodes = ifelse(direction == "Score_COS_node1",
                          paste(node1, node2, sep = "/"),
                          paste(node2, node1, sep = "/"))
) #ifelse() Verifica el valor de la columna direction. Si es "Score_COS_node1", 
#concatena node1 y node2. Si es "Score_COS_node2", concatena node2 y node1.

combined_interactions <- combined_interactions %>%
  separate(combined_nodes, into = c("nodeA", "nodeB"), sep = "/")# Ahora nodo A toma la posición del nodo de interes

# Eliminar filas donde nodeA y nodeB sean iguales, es decir, singuletes
combined_interactions <- combined_interactions[combined_interactions$nodeA != combined_interactions$nodeB, ]

score_COS <- combined_interactions %>%
  select(nodeA, cosine_score) # Se filtra el valor Coseno según el nodo de interes (Unido a un nodo anotado)

# Merge dataframes by id
BASE <- merge(score_MSQ, score_COS, by.x = "id", by.y = "nodeA", all = TRUE) # Note that the MQS value is repeated when there is more than one COS, which is exactly what was sought


### Retrieving the third parameter, number of unannotated nodes
score_no_MSQ <- df_nodes %>%
  select(id, MQScore, UniqueFileSources) %>%
  filter(is.na(MQScore) | MQScore == "")%>%
  select(-MQScore)

df_nodes_unannotated <- score_no_MSQ %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>% #ACÁ TAMBIÉN HAY QUE CAMBIAR EL IDENTIFACOR DE MUESTRA
  group_by(UniqueFileSources)

# Count id per sample
conteo <- table(df_nodes_unannotated$UniqueFileSources)
df_conteo <- as.data.frame(conteo)
total_freq <- sum(df_conteo$Freq)

# Calcular el porcentaje de cada frecuencia
df_conteo$Percentage <- (df_conteo$Freq / total_freq)*100 # se mantiene en porcentaje por la escala con los otros parámetros

df_conteo <- df_conteo %>% select(-Freq)

######Shannon index for MQS y Cos

# Crear la carpeta CSV_IQ si no existe
if (!dir.exists("CSV_IQ")) {
  dir.create("CSV_IQ")
}

# Merge dataframes by id
BASE <- merge(score_MSQ, score_COS, by.x = "id", by.y = "nodeA", all = TRUE)
# Note that the MQS value is repeated when there is more than one COS, which is exactly what was sought

# Split samples and extend the dataframe
df_separated <- BASE %>%
  separate_rows(UniqueFileSources, sep = "\\|") %>%
  group_by(UniqueFileSources)

# Define the strata
strata_MSQ <- list(
  "0.01-0.05" = c(0.01, 0.05),
  "0.05-0.10" = c(0.05, 0.10),
  "0.10-0.15" = c(0.10, 0.15),
  "0.15-0.20" = c(0.15, 0.20),
  "0.20-0.25" = c(0.20, 0.25),
  "0.25-0.30" = c(0.25, 0.30),
  "0.30-0.35" = c(0.30, 0.35),
  "0.35-0.40" = c(0.35, 0.40),
  "0.40-0.45" = c(0.40, 0.45),
  "0.45-0.50" = c(0.45, 0.50),
  "0.50-0.55" = c(0.50, 0.55),
  "0.55-0.60" = c(0.55, 0.60),
  "0.60-0.65" = c(0.60, 0.65),
  "0.65-0.70" = c(0.65, 0.70),
  "0.70-0.75" = c(0.70, 0.75),
  "0.75-0.80" = c(0.75, 0.80),
  "0.80-0.85" = c(0.80, 0.85),
  "0.85-0.90" = c(0.85, 0.90),
  "0.90-0.95" = c(0.90, 0.95),
  "0.95-1.00" = c(0.95, 1.00)
)

strata_COS <- list(
  "0.70-0.75" = c(0.70, 0.75),
  "0.75-0.80" = c(0.75, 0.80),
  "0.80-0.85" = c(0.80, 0.85),
  "0.85-0.90" = c(0.85, 0.90),
  "0.90-0.95" = c(0.90, 0.95),
  "0.95-1.00" = c(0.95, 1.00)
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
  MQScore_counts <- stratify_and_count_unique(individual_df, "MQScore", strata_MSQ)
  
  # Count values in each strata for cosine_score
  cosine_score_counts <- stratify_and_count(individual_df, "cosine_score", strata_COS)
  
  # Calculate Shannon Index for MQScore
  MQScore_shannon <- calculate_shannon_index(MQScore_counts)
  
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


################

# Define the new number of samples per plot
samples_per_plot <- 10   # Cambiar este valor según sea necesario

# Recalculate the number of plots needed based on the new value
num_plots <- ceiling(nrow(IQ_Summary) / samples_per_plot)

# Crear la carpeta IQ_Plot si no existe
if (!dir.exists("IQ_Plot")) {
  dir.create("IQ_Plot")
}

# Generar los gráficos 
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

# resumen de las 10 muestras con los dos menores índices Shannon
IQ_Summary_top <- IQ_Summary %>%
  mutate(IQ= Hmsq + Hcos - Nsna) %>%
  arrange(IQ) %>%
  head(10)

# Crear el gráfico resumen con los parámetros
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

# Guardar el gráfico resumen en la carpeta IQ_Plot
ggsave("IQ_Plot/Parameters_summary.png", plot = p_summary, width = 15, height = 8)


# Crear el gráfico resumen de IQ
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

# Guardar el gráfico resumen en la carpeta IQ_Plot
ggsave("IQ_Plot/IQ_Plot_summary.png", plot = p_summary, width = 15, height = 8)

###########################################################################


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

# Resultado BiNi
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

# Graficamos
IQ_BiNI_Plot<-ggplot(df_combined, aes(x = solvent, y = IQ, color = solvent)) +
  geom_point(size = 3) +
  facet_wrap(~ Strain, scales = "free_y") +
  geom_hline(data = biNI_values, aes(yintercept = BiNI_value), linetype = "dashed", color = "red") +
  labs(title = "IQ por Strain y Solvent con Línea de BiNI",
       x = "Solvent",
       y = "IQ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ # Rota etiquetas del eje x si es necesario
  ylim(1, 4) # Establece los límites del eje y 

# Guardar el gráfico resumen en la carpeta IQ_Plot
ggsave("IQ_Plot/IQ_BiNI_Plot.png", plot = IQ_BiNI_Plot, width = 15, height = 8)
# Save the updated Shannon Index summary dataframe to a CSV file in the CSV_IQ folder
write_xlsx(df_combined, path = "CSV_IQ/2_IQ_BINI.xlsx")









###########################################################################
#1- sE DEBE TOMAR UN TOP N DE LOS NO ANOTADOS QUE TENGA EL MENOR MSQ Y EL MENOR VALOR ASOCIADO CON LA FAMILIA A LA QUE PERTENECE
#2- A ESOS FEATURES, DETERMINAR SIRIUS, OBTENER UNA TABLA DE PUNTAJES DE SIRIUS Y CSI E IMPRIMIR UN CSV