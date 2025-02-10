# ---- calculating network properties ----
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(vegan)

## ---- functions ----
build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')
  
  # Step 3: Create the matrix with specific row and column species
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows
  
  # Initialize an empty matrix
  interaction_matrix <- matrix(0, nrow = length(species_to), ncol = length(species_from),
                               dimnames = list(species_to, species_from))
  
  # Populate the matrix with aggregated weights
  for (i in 1:nrow(aggregated_data)) {
    row <- aggregated_data$node_to[i]    # Rows represent 'node_to' species
    col <- aggregated_data$node_from[i]  # Columns represent 'node_from' species
    interaction_matrix[row, col] <- aggregated_data$weight[i]
  }
  
  return(interaction_matrix)
}
## ---- parameters ----
emln_id <- 60

## ---- run ----
### ---- load matrices ----
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict
for (net in emln_id) {
  
  # Load matrices
  d <- load_emln(net)  # Use 'net' instead of 'emln_id'
  graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
  A_l <- d$extended
  num_layers <- length(graph_list)  # Define before the loop
  
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
      
      ## ---- calculating bray-curtis ----
      
      # Get all unique row and column names
      all_rows <- union(rownames(A), rownames(P))
      all_cols <- union(colnames(A), colnames(P))
      
      # Expand both matrices to include all rows and columns, filling missing values with 0
      A_expanded <- matrix(0, nrow = length(all_rows), ncol = length(all_cols), 
                           dimnames = list(all_rows, all_cols))
      P_expanded <- matrix(0, nrow = length(all_rows), ncol = length(all_cols), 
                           dimnames = list(all_rows, all_cols))
      
      # Fill the matrices with existing values
      A_expanded[rownames(A), colnames(A)] <- A
      P_expanded[rownames(P), colnames(P)] <- P
      
      # Convert matrices to vectors for comparison
      A_vector <- as.vector(A_expanded)
      P_vector <- as.vector(P_expanded)
      
      # Compute Bray-Curtis similarity (1 - dissimilarity)
      bray_curtis_similarity <- 1 - vegdist(rbind(A_vector, P_vector), method = "bray")
      
      # Print the similarity score
      print(bray_curtis_similarity)
      
      # Add these values to the results table
      results <- rbind(results, data.frame(emln_id = emln_id,
                                           train_layer = layers_to_train,
                                           test_layer = layer_to_predict,
                                           bray_curtis_similarity = bray_curtis_similarity))  
      
      
      
    }}}
