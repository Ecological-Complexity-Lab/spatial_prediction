# ---- calculating network dimensions ----
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())
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
# Load matrices
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
A_l <- d$extended

# Total number of layers
num_layers <- length(graph_list)

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

# View results
View(results)

results %>% write_csv('result_netsize_canaries_site_scale.csv')

# add to main results
result_summary <- read_csv('working_df_all_itr_60_binary_names.csv')
result_summary <- read_csv('working_df_site_scaled_evaluators_distance.csv') # scaled

result_summary <- result_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

result_summary %>% write_csv('result_netsize_canaries_distance_names_site_scaled.csv')

## ---- checking correlations ----
### ---- 1 off-diagonal no diagonal ----
# use half the matrix ('cause 1 <- 2 same as 2 <- 1)

canary_results_1off_no_diag <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer)

# if we want to consider all data points, replace canary_results_1off_no_diag with result_summary hereafter
df_long_1off <- canary_results_1off_no_diag %>%
  pivot_longer(
    cols = c(size_P, density_P, size_C, density_C),
    names_to = "measure_type",
    values_to = "measure_value"
  )

# For each variable, compute correlation with f1_score:
cor_table <- df_long_1off %>%
  group_by(measure_type) %>%
  summarise(
    cor_value = cor(f1_score, measure_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, measure_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table

cor_table_annot <- cor_table %>%
  mutate(
    # round correlation to 3 decimals, no scientific notation
    r_fmt  = formatC(cor_value, format = "f", digits = 2),
    # round p-value to 4 decimals, no scientific notation
    p_fmt  = formatC(p_value,  format = "f", digits = 4),
    label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
  )

# Create a named vector for renaming facets
facet_labels <- c(
  "size_C" = "Size of matrix C",
  "density_C" = "Density of matrix C",
  "size_P" = "Size of matrix P",
  "density_P" = "Density of matrix P"
)

ggplot(df_long_1off, aes(x = measure_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +
  facet_wrap(
    ~ measure_type,
    scales   = "free_x",
    labeller = as_labeller(facet_labels)  # Use the named vector
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  geom_text(
    data    = cor_table_annot,
    aes(label = label_text),
    x       = Inf,
    y       = Inf,
    hjust   = 1.1,
    vjust   = 1.2,
    size    = 3.2,
    color   = "black"
  ) +
  labs(
    x = "Network feature",
    y = "F1 score",
    title = "F1 vs. network measures - 1 off-diagonal"
  ) +
  theme_minimal() +
  tme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )
