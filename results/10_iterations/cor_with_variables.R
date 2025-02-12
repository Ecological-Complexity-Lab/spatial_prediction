# ---- calculating network properties ----
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(vegan)
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

## ---- run ----
### ---- load matrices ----
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict. this loop gives us the bray-curtis similarity between the interaction composition
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


# View final results
print(results)

 ## ---- plot bc ----
result %>% write_csv('result_bray_curtis_canaries.csv')

result_summary <- read_csv('result_summary_canary_with_distance.csv')

result_summary <- result_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table. bray_curtis is for total similarity

# use half the matrix

canary_results_filtered <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

# Compute correlation
correlation <- cor.test(canary_results_filtered$f1_score, canary_results_filtered$bray_curtis_similarity, 
                        use = "complete.obs", method = "pearson")

# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # Use scientific notation if needed
label_text <- paste0("r = ", r_value, ", p = ", p_value)

ggplot(canary_results_filtered, aes(x = bray_curtis_similarity, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +  # Trendline
  labs(x = "Overall Bray-Curtis similarity",
       y = "F1 score",
       title = "F1 score vs. overall bray_curtis similarity for Canary results") +
  tme +  # Corrected theme
  annotate("text",
           x = max(canary_results_filtered$bray_curtis_similarity) * 0.9,  # Position text dynamically
           y = max(canary_results_filtered$f1_score) * 1,  
           label = label_text,
           size = 3,
           color = "black")


## ---- Jaccard for plants, pollinators and links ----
results_jaccard <- data.frame()

for (net in emln_id) {
  
  # Load matrices
  d <- load_emln(net)  # Use 'net' instead of 'emln_id'
  graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
  A_l <- d$extended
  num_layers <- length(graph_list)  # Define before the loop
  
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      
      A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
      P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
      
      # 1) Jaccard pollinators
      poll_train <- rownames(A)[ rowSums(A) > 0 ]
      poll_test  <- rownames(P)[ rowSums(P) > 0 ]
      intersection_poll <- length(intersect(poll_train, poll_test))
      union_poll        <- length(union(poll_train, poll_test))
      jaccard_poll <- if (union_poll == 0) NA else intersection_poll / union_poll
      
      # 2) Jaccard plants
      plants_train <- colnames(A)[ colSums(A) > 0 ]
      plants_test  <- colnames(P)[ colSums(P) > 0 ]
      intersection_plants <- length(intersect(plants_train, plants_test))
      union_plants        <- length(union(plants_train, plants_test))
      jaccard_plants <- if (union_plants == 0) NA else intersection_plants / union_plants
      
      # 3) Jaccard edges
      pairs_train <- which(A > 0, arr.ind = TRUE)
      pairs_train_strings <- apply(pairs_train, 1, function(rc) {
        paste(rownames(A)[rc[1]], colnames(A)[rc[2]], sep = "_")
      })
      
      pairs_test <- which(P > 0, arr.ind = TRUE)
      pairs_test_strings <- apply(pairs_test, 1, function(rc) {
        paste(rownames(P)[rc[1]], colnames(P)[rc[2]], sep = "_")
      })
      intersection_edges <- length(intersect(pairs_train_strings, pairs_test_strings))
      union_edges        <- length(union(pairs_train_strings, pairs_test_strings))
      jaccard_edges <- if (union_edges == 0) NA else intersection_edges / union_edges
      
      # store the result
      results_jaccard <- rbind(
        results_jaccard,
        data.frame(
          emln_id = net,
          train_layer = layers_to_train,
          test_layer  = layer_to_predict,
          jaccard_pollinators = jaccard_poll,
          jaccard_plants      = jaccard_plants,
          jaccard_edges       = jaccard_edges
        )
      )
    }
  }
}

head(results_jaccard)

## ---- analyse ----

results_jaccard %>% write_csv('result_jaccard_canaries.csv')

result_summary <- read_csv('result_summary_canary_with_distance.csv')

result_summary <- result_summary %>%
  left_join(results_jaccard, by = c("train_layer", "test_layer")) # add to results table. bray_curtis is for total similarity

# use half the matrix

canary_results_filtered <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

# Example data frame: canary_results_filtered
# with columns: f1_score, jaccard_pollinators, jaccard_plants, jaccard_edges, etc.

df_long <- canary_results_filtered %>%
  pivot_longer(
    cols = c(jaccard_pollinators, jaccard_plants, jaccard_edges),
    names_to = "jaccard_type",
    values_to = "jaccard_value"
  )

ggplot(df_long, aes(x = jaccard_value, y = f1_score, color = jaccard_type)) +
  geom_point(alpha = 0.6, size = 2) +   # scatter points
  geom_smooth(method = "lm", se = FALSE) +  # linear fit per jaccard_type
  labs(
    x = "Jaccard Similarity",
    y = "F1 Score",
    title = "F1 vs. Jaccard Measures"
  ) +
  theme_minimal()

# only for off-diagonals

canary_results_offs <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer > test_layer)

df_long_offs <- canary_results_offs %>%
  pivot_longer(
    cols = c(jaccard_pollinators, jaccard_plants, jaccard_edges),
    names_to = "jaccard_type",
    values_to = "jaccard_value"
  )

ggplot(df_long_offs, aes(x = jaccard_value, y = f1_score, color = jaccard_type)) +
  geom_point(alpha = 0.6, size = 2) +   # scatter points
  geom_smooth(method = "lm", se = FALSE) +  # linear fit per jaccard_type
  labs(
    x = "Jaccard Similarity",
    y = "F1 Score",
    title = "F1 vs. Jaccard Measures"
  ) +
  theme_minimal()

ggplot(df_long_offs, aes(x = jaccard_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  labs(x = "Jaccard similarity", y = "F1 score") +
  facet_wrap(~ jaccard_type, scales = "free_x") +
  theme_minimal() +
  tme

### ---- calculating correlation ----
# For each jaccard_type, compute correlation with f1_score:
cor_table <- df_long_offs %>%
  group_by(jaccard_type) %>%
  summarise(
    cor_value = cor(f1_score, jaccard_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, jaccard_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table


# Convert cor_table to a format suitable for annotation:
cor_table_annot <- cor_table %>%
  mutate(
    label_text = paste0("r = ", round(cor_value, 3),
                        "\np = ", formatC(p_value, digits = 3))
  )

ggplot(df_long_offs, aes(x = jaccard_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  labs(x = "Jaccard similarity", y = "F1 score") +
  facet_wrap(~ jaccard_type, scales = "free_x") +
  theme_minimal() +
  # Add correlation label in each facet
  geom_text(
    data = cor_table_annot,
    aes(label = label_text), 
    x = Inf, y = Inf,  # place near top-right corner
    hjust = 1.1, vjust = 1.2,
    color = "black",
    size = 3
  ) + tme

