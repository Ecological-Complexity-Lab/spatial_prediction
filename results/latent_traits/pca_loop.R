# say your original data is in a matrix or data.frame `X`
# rows = species, columns = ecological variables

## ---- run ----
# run the algorithm with 1 iteration
# ---- softImpute for predicting removed links in empirical networks - all layer combinations ----
# including removal of zeros and ones. no bootstrapping yet
# in this code we remove zeros and ones and get the predicted values (not labels but values based on the svd from the softImpute fit function) in a dataframe for further exploration. We get results for a range of k and lambdas.
## ---- to do ----
#   CHECK IF COORDINATES IMPLY WHAT SHOULD BE INSTEAD OF NA - AND THEN WHAT IS THE MEANING OF REMOVING LINKS?
# check that there are no duplicates
# bootstrapping
# add aggregation (same format, different file)
# fix combined_results to include all k and lambda values for non-removed links
# ---- results so far -----
# results for all links and layer to layer predictions in combined_results_0.2_rem_values_nonbinary_all_edges2.csv.
## ---- load libraries ----
library(softImpute)
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(factoextra)
library(tidyverse)

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
layers_to_train <- 1
layer_to_predict <- 7
prop_ones_to_remove <- 0.2
# prop_zeros_to_remove <- 0.2
n_sim <- 1
is_binary <- 0
set.seed(42)

## ---- run ----
### ---- load matrices ----

# Load matrices
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
A_l <- d$extended

# aggregate to island scale
# Extract numeric layer numbers
A_l <- A_l %>%
  mutate(layer_num = as.numeric(gsub("layer_", "", layer_from))) %>%
  mutate(aggregated_layer = ifelse(layer_num %% 2 == 1, 
                                   paste0("layer_", layer_num, "_", layer_num + 1),
                                   paste0("layer_", layer_num - 1, "_", layer_num)))3
# Total number of layers
num_layers <- length(unique(A_l$layer_from))

# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layer_to_predict in 1:num_layers) {
  # Build the layer to predict matrix P
  P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
  
### ---- pca ----
  svd_res <- svd(P)
  
  # svd_res$u is n_poll × n_poll
  # svd_res$v is n_plant × n_plant
  # svd_res$d is length = min(n_poll,n_plant)
  
  # Compute “scores” = U %*% Σ  and  V %*% Σ
  poll_scores  <- svd_res$u %*% diag(svd_res$d)    # rows = pollinators
  plant_scores <- svd_res$v %*% diag(svd_res$d)    # rows = plants
  
  poll_df <- as_tibble(
    poll_scores[,1:2],
    .name_repair = ~ c("PC1","PC2")
  ) %>%
    add_column(species = rownames(P_original), .before=1)
  
  plant_df <- as_tibble(
    plant_scores[,1:2],
    .name_repair = ~ c("PC1","PC2")
  ) %>%
    add_column(species = colnames(P_original), .before=1)
  
  # Pollinators
  ggplot(poll_df, aes(x = PC1, y = PC2, label = species)) +
    geom_point() +
    geom_text_repel(size = 3) +
    labs(title = "SVD‐ordination: pollinators",
         x     = "Axis 1",
         y     = "Axis 2") +
    theme_minimal() + tme
  
  ggplot(poll_df, aes(x = PC1, y = PC2, label = species)) +
    geom_point() +
    geom_text_repel(
      size         = 3,
      max.overlaps = Inf,      # ← allow *all* labels
      box.padding  = 0.5,      # tweak these if labels still clash
      point.padding= 0.3,
      force        = 1.0       # default is 1; increase to push labels further apart
    ) +
    labs(
      title = "SVD‑ordination: pollinators",
      x     = "Axis 1",
      y     = "Axis 2"
    ) +
    theme_minimal() + tme
  
  # Plants
  ggplot(plant_df, aes(x = PC1, y = PC2, label = species)) +
    geom_point() +
    geom_text_repel(size = 4) +
    labs(title = "SVD‐ordination: plants",
         x     = "Axis 1",
         y     = "Axis 2") +
    theme_minimal() + tme
  
  ggplot(plant_df, aes(x = PC1, y = PC2, label = species)) +
    geom_point() +
    geom_text_repel(
      size         = 3,
      max.overlaps = Inf,      # ← allow *all* labels
      box.padding  = 0.5,      # tweak these if labels still clash
      point.padding= 0.3,
      force        = 1.0       # default is 1; increase to push labels further apart
    ) +
    labs(
      title = "SVD‑ordination: plants",
      x     = "Axis 1",
      y     = "Axis 2"
    ) +
    theme_minimal() + tme
}