# ---- calculating network dimensions ----
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)

## ---- parameters ----
emln_id <- 60

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
      
      node_to <- rownames(P) # for the results
      node_from <- colnames(P)
      
      ### ---- creating a combined matrix C ----
      all_row_ids <- unique(c(rownames(A), rownames(P)))
      all_col_ids <- unique(c(colnames(A), colnames(P)))
      C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                  dimnames = list(all_row_ids, all_col_ids))
      
      # Place A into C
      C[rownames(A), colnames(A)] <- A
      
      # Place P into C
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
      
      # Compute matrix properties
      nrow_A <- nrow(A)
      nrow_P <- nrow(P)
      nrow_C <- nrow(C)
      ncol_A <- ncol(A)
      ncol_P <- ncol(P)
      ncol_C <- ncol(C)
      size_A <- length(A)
      size_P <- length(P)
      size_C <- length(C)
      
      # Calculate density for A, P, and C
      density_A <- sum(A > 0) / length(A)
      density_P <- sum(P > 0) / length(P)
      density_C <- sum(C > 0) / length(C)
      
      # Add these values to the results table
      results <- rbind(results, data.frame(emln_id = emln_id,
                                           train_layer = layers_to_train,
                                           test_layer = layer_to_predict,
                                           nrow_A = nrow_A,
                                           ncol_A = ncol_A,
                                           size_A = size_A,
                                           density_A = density_A,   # Added density
                                           nrow_P = nrow_P,
                                           ncol_P = ncol_P,
                                           size_P = size_P,
                                           density_P = density_P,   # Added density
                                           nrow_C = nrow_C,
                                           ncol_C = ncol_C,
                                           size_C = size_C,
                                           density_C = density_C))  # Added density
      
    
    }
  }
}

# add this
data(olesen2002flores) # Check out the help for information on this data set!
olesen2002flores_binary <- 1 * (olesen2002flores > 0) # Make the data binary (unweighted)
I <- nrow(olesen2002flores_binary) # Number of lower level species (e.g., hosts, plants)
J <- ncol(olesen2002flores_binary) # Number of higher level species (e.g., parasites, pollinators)
S <- I + J # Total number of species, aka: Network size
L <- sum(olesen2002flores_binary > 0) # Number of edges in the network
A_i <- rowSums(olesen2002flores_binary) # The degree of hosts
A_j <- colSums(olesen2002flores_binary) # The degree of parasites
C <- L / (I * J) # Connectance
# Clustering coefficient higher level (the number of realized links 
# divided by the number of possible links for each species)
cc_high <- colSums(olesen2002flores_binary) / nrow(olesen2002flores_binary) 
# Clustering coefficient lower level 
# (the number of realized links divided by the number of possible links for each species)
cc_low <- rowSums(olesen2002flores_binary) / ncol(olesen2002flores_binary) 

# View results
View(results)

# Select only the relevant columns from result_summary
result_summary_subset <- result_summary[, c("emln_id", "train_layer", "test_layer", "f1_score")]

# Merge only the f1_score_sigm column into results
results <- merge(results, result_summary_subset, 
                 by = c("emln_id", "train_layer", "test_layer"), 
                 all.x = TRUE)

# View the updated results table
View(results)

# Filter results for emln_id = 60
#canary_results <- subset(results, emln_id == 60)
canary_results <- results
## ---- checking correlations ----
# Compute correlation
correlation <- cor.test(canary_results$f1_score, canary_results$size_C, 
                        use = "complete.obs", method = "pearson")

# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # Use scientific notation if needed
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Create scatter plot
ggplot(canary_results, aes(x = size_C, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +  # Trendline
  labs(x = "Size of C",
       y = "F1 score",
       title = "F1 score vs. size of C for Canary results") +
  tme +  # Corrected theme
  annotate("text",
           x = max(canary_results$size_C) * 0.9,  # Position text dynamically
           y = max(canary_results$f1_score) * 1,  
           label = label_text,
           size = 3,
           color = "black")

## ---- now with one off-diagonal ----

canary_results_filtered <- canary_results %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

# Compute correlation
correlation <- cor.test(canary_results_filtered$f1_score, canary_results_filtered$density_C, 
                        use = "complete.obs", method = "pearson")

# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # Use scientific notation if needed
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Create scatter plot
ggplot(canary_results_filtered, aes(x = density_C, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +  # Trendline
  labs(x = "Density of C",
       y = "F1 score",
       title = "F1 score vs. density of C for Canary results (1 diagonal)") +
  tme +  # Corrected theme
  annotate("text",
           x = max(canary_results_filtered$density_C) * 0.9,  # Position text dynamically
           y = max(canary_results_filtered$f1_score) * 1,  
           label = label_text,
           size = 3,
           color = "black")
