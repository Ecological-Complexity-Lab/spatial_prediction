# ---- Jaccard island scale ----
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

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

emln_id <- 60
## ---- Jaccard for plants, pollinators and links ----

### ---- load matrices ----
# Initialize a data frame to store combined results for all layer combinations
results_jaccard <- data.frame()

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
                                   paste0("layer_", layer_num - 1, "_", layer_num)))

# Aggregate data
aggregated_df <- A_l %>%
  group_by(aggregated_layer, node_from, node_to, type) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  mutate(layer_from = aggregated_layer, layer_to = aggregated_layer) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# Generate new layer names
unique_layers <- unique(aggregated_df$layer_from)  # Get unique aggregated layer names
new_layer_names <- paste0("layer_", seq_along(unique_layers))  # Generate new names (layer_1, layer_2, ...)

# Create a mapping table
layer_mapping <- data.frame(original_layer = unique_layers, new_layer = new_layer_names)

# Apply renaming in aggregated_df
aggregated_df <- aggregated_df %>%
  left_join(layer_mapping, by = c("layer_from" = "original_layer")) %>%
  mutate(layer_from = new_layer, layer_to = new_layer) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# View updated aggregated_df
print(aggregated_df)

# Total number of layers
num_layers <- length(unique(aggregated_df$layer_from))

# if we want Jaccard for site scale:
#num_layers <- length(unique(A_l$layer_from))
  
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
          emln_id = emln_id,
          train_layer = layers_to_train,
          test_layer  = layer_to_predict,
          jaccard_pollinators = jaccard_poll,
          jaccard_plants      = jaccard_plants,
          jaccard_edges       = jaccard_edges
        )
      )
    }
  }


head(results_jaccard)

#results_jaccard %>% write_csv('result_jaccard_canaries_island_scaled.csv')

results_jaccard <- read.csv('result_jaccard_canaries_island_scale.csv')
results_jaccard <- read.csv('result_jaccard_canaries_site_scaled.csv')

result_summary <- read_csv('working_df_island_distance_fidelity.csv')
result_summary <- read_csv('working_df_all_itr_60_binary_scaled_site_names.csv') # site scale, centered version
result_summary <- read.csv('working_df_islands_scaled_evaluators_distance.csv') # island scale, centered

result_summary <- result_summary %>%
  left_join(results_jaccard, by = c("train_layer", "test_layer")) # add to results table
# 
# result_summary <- read_csv('result_netsize_canaries_distance_names_site_scaled.csv') # site scale, centered version
# result_summary <- result_summary %>%
#   left_join(results_jaccard, by = c("train_layer", "test_layer")) # add to results table
# result_summary %>% write_csv('result_netsize_canaries_distance_names_jaccard_site_scaled.csv')

# result_summary <- read_csv('working_df_islands_scaled_evaluators_distance.csv') # island scale, centered version
# result_summary <- result_summary %>%
#   left_join(results_jaccard, by = c("train_layer", "test_layer")) # add to results table
# result_summary %>% write_csv('result_canaries_distance_names_jaccard_island_scaled.csv')
## ---- plot ----
### ---- only 1 off-diagonal and diagonal ----
# use half the matrix ('cause 1 <- 2 same as 2 <- 1)

canary_results_1off <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

df_long_1off <- canary_results_1off %>%
  pivot_longer(
    cols = c(jaccard_pollinators, jaccard_plants, jaccard_edges),
    names_to = "jaccard_type",
    values_to = "jaccard_value"
  )

# For each jaccard_type, compute correlation with f1_score:
cor_table <- df_long_1off %>%
  group_by(jaccard_type) %>%
  summarise(
    cor_value = cor(f1_score, jaccard_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, jaccard_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table

cor_table_annot <- cor_table %>%
  mutate(
    # round correlation to 3 decimals, no scientific notation
    r_fmt  = formatC(cor_value, format = "f", digits = 2),
    # round p-value to 4 decimals, no scientific notation
    p_fmt  = formatC(p_value,  format = "f", digits = 2),
    label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
  )

ggplot(df_long_1off, aes(x = jaccard_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  facet_wrap(
    ~ jaccard_type,
    scales   = "free_x",             # or "free" if you want x & y free
    labeller = as_labeller(type_labels)  # rename facets
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), ) +
  # Annotation: place correlation in top-right corner of each facet
  geom_text(
    data    = cor_table_annot,
    aes(label = label_text),
    x       = Inf,
    y       = Inf,
    hjust   = 1.1,  # move left from right edge
    vjust   = 2.2,  # move down from top edge
    size    = 3.2,
    color   = "black"
  ) +
  labs(
    x = "Jaccard similarity",
    y = "F1 score",
    title = "F1 vs. Jaccard measures - 1 off-diagonal + diag"
  ) +
  theme_minimal() +
  tme +
  theme(
    # Add a black frame around facet labels with thickness
    #strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )

### ---- 1 off-diagonal no diagonal ----
# use half the matrix ('cause 1 <- 2 same as 2 <- 1)

canary_results_1off_no_diag <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer)

df_long_1off <- canary_results_1off_no_diag %>%
  pivot_longer(
    cols = c(jaccard_pollinators, jaccard_plants, jaccard_edges),
    names_to = "jaccard_type",
    values_to = "jaccard_value"
  )

# For each jaccard_type, compute correlation with f1_score:
cor_table <- df_long_1off %>%
  group_by(jaccard_type) %>%
  summarise(
    cor_value = cor(f1_score, jaccard_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, jaccard_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table

cor_table_annot <- cor_table %>%
  mutate(
    # round correlation to 3 decimals, no scientific notation
    r_fmt  = formatC(cor_value, format = "f", digits = 2),
    # round p-value to 4 decimals, no scientific notation
    p_fmt  = formatC(p_value,  format = "f", digits = 2),
    label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
  )

ggplot(df_long_1off, aes(x = jaccard_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  facet_wrap(
    ~ jaccard_type,
    scales   = "free_x"
    #,             # or "free" if you want x & y free
    #labeller = as_labeller(type_labels)  # rename facets
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), ) +
  # Annotation: place correlation in top-right corner of each facet
  geom_text(
    data    = cor_table_annot,
    aes(label = label_text),
    x       = Inf,
    y       = Inf,
    hjust   = 1.1,  # move left from right edge
    vjust   = 1.2,  # move down from top edge
    size    = 3.2,
    color   = "black"
  ) +
  labs(
    x = "Jaccard similarity",
    y = "F1 score",
    title = "F1 vs. Jaccard measures - 1 off-diagonal"
  ) +
  theme_minimal() +
  tme +
  theme(
    # Add a black frame around facet labels with thickness
    #strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )

### ---- 2 off-diagonals no diagonal ----

canary_results_diags <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer != test_layer)

df_long_1off <- canary_results_diags %>%
  pivot_longer(
    cols = c(jaccard_pollinators, jaccard_plants, jaccard_edges),
    names_to = "jaccard_type",
    values_to = "jaccard_value"
  )

# For each jaccard_type, compute correlation with f1_score:
cor_table <- df_long_1off %>%
  group_by(jaccard_type) %>%
  summarise(
    cor_value = cor(f1_score, jaccard_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, jaccard_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table

cor_table_annot <- cor_table %>%
  mutate(
    # round correlation to 3 decimals, no scientific notation
    r_fmt  = formatC(cor_value, format = "f", digits = 2),
    # round p-value to 4 decimals, no scientific notation
    p_fmt  = formatC(p_value,  format = "f", digits = 2),
    label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
  )

ggplot(df_long_1off, aes(x = jaccard_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  facet_wrap(
    ~ jaccard_type,
    scales   = "free_x"
    #,             # or "free" if you want x & y free
    #labeller = as_labeller(type_labels)  # rename facets
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), ) +
  # Annotation: place correlation in top-right corner of each facet
  geom_text(
    data    = cor_table_annot,
    aes(label = label_text),
    x       = Inf,
    y       = Inf,
    hjust   = 1.1,  # move left from right edge
    vjust   = 1.2,  # move down from top edge
    size    = 3.2,
    color   = "black"
  ) +
  labs(
    x = "Jaccard similarity",
    y = "F1 score",
    title = "F1 vs. Jaccard measures - 1 off-diagonal"
  ) +
  theme_minimal() +
  tme +
  theme(
    # Add a black frame around facet labels with thickness
    #strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )

