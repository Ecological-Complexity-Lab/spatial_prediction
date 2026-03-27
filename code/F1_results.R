# 2 figures with F1 - separate script
# here we obtain a confusion matrix and figures based on the previous, 
# less restrictive evaluation metod used before revising the evaluation to be based on the F0.5 score.

# Load necessary libraries
library(ggplot2)

# Load the data
original_run_combined_results <- readRDS(results_file) # produced in the main script


# fig. 2c - heatmap red (island_heatmap_f1.pdf) ----


# after reading or producing the results, filter these:
or_combined_results <- original_run_combined_results %>% 
  filter(k == 2) %>% 
  filter(!(input_lambda %in%  c(1, 5, 50, 100)))

# convert negatives to zeros
or_df <- or_combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

or_df_removed <- or_df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

## calc theashold ----

or_df_avg <- or_df %>% # 1) filter & prep
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  ) %>% # 2) expand to one row per threshold
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy= (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) /
      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)#,
    #.groups = "drop"
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup() %>% # 3) average across emln_id/layer combos and pivot long
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall,
      f1_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold,
               names_to  = "metric",
               values_to = "value")

or_df_wide <- or_df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# # 2) find the threshold with the optimal f1 score
# or_best_discrete <- or_df_wide %>%
#   slice_max(f1_score, n = 1)

# find threshold based on f1-balanced accuracy trade-off
or_best_discrete <- or_df_wide

or_best_discrete <- or_df_wide %>%
  mutate(absdiff = abs(f1_score - balanced_accuracy)) %>%
  slice_min(absdiff, n = 1)

# results:
or_best_discrete_threshold <- or_best_discrete$threshold
or_best_discrete_threshold


or_df_removed <- or_df %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > or_best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

# calculate summary statistics
or_result_summary <- or_df_removed %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
      sum((original_links   - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
    nnse = 1 / (2 - nse)
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()

or_result_summary <- or_result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Single location",
    train_layer != test_layer ~ "Added location"
  ))


or_results <- data.frame()

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    ### creating a combined matrix C 
    all_row_ids <- unique(c(rownames(A), rownames(P)))
    all_col_ids <- unique(c(colnames(A), colnames(P)))
    C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                dimnames = list(all_row_ids, all_col_ids))
    
    # Place A into C
    C[rownames(A), colnames(A)] <- A
    
    # Place P into C
    if (layers_to_train != layer_to_predict){
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    } else {
      # if this predicts using the same layer, don't sum it to itself
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            (C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])/2)
    }
    
    
    # make them all binary for count
    A[A>0] <- 1
    P[P>0] <- 1
    C[C>0] <- 1
    
    # Compute matrix properties
    nrow_A <- nrow(A)
    nrow_P <- nrow(P)
    nrow_C <- nrow(C)
    ncol_A <- ncol(A)
    ncol_P <- ncol(P)
    ncol_C <- ncol(C)
    size_A <- nrow(A) + ncol(A)
    size_P <- nrow(P) + ncol(P)
    size_C <- nrow(C) + ncol(C)
    
    # Calculate density for A, P, and C
    density_A <- sum(A > 0) / length(A)
    density_P <- sum(P > 0) / length(P)
    density_C <- sum(C > 0) / length(C)
    
    # Add these values to the results table
    or_results <- rbind(or_results, data.frame(emln_id = emln_id,
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

# view results
summary(or_results)


# add them to the results table
or_result_summary <- or_result_summary %>%
  left_join(or_results, by = c("train_layer", "test_layer")) # add to results table

# summerize (table ST1)

# Add island names to main table
or_result_summary <- or_result_summary %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)

# now we can summarise
or_df_summary <- or_result_summary %>%
  group_by(test_layer_name) %>%
  summarise(
    size_P   = mean(size_P, na.rm = TRUE),
    density_P  = mean(density_P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(size_P))  # change to asc() for smallest first


or_result_summary_island <- or_result_summary

# Add to main table
or_result_summary_island <- or_result_summary_island %>%
  left_join(
    distance_island_table,
    by = c("train_layer_name" = "from", "test_layer_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_layer_name == test_layer_name,
                               0,              # distance = 0 if same site
                               avg_distance_km))   # otherwise, keep joined distance

or_island_heatmap_f1 <- 
  ggplot(or_result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = or_result_summary_island[or_result_summary_island$train_layer == or_result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added location", y = "Predicted location", fill = "F0.5 score") +
  theme_minimal() +
  theme(
    text = element_text(size = 9),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(or_island_heatmap_f1)

pdf(
  file   = "results/paper_figs/island_heatmap_f1_as_cutoff.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(or_island_heatmap_f1)
dev.off()     # close the file

# the figures:
# parallel to fig. 3a - map_missing_links.pdf ----
# calculate best threshold


or_df <- or_df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

or_df_island <- or_df %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )

table(or_df_island$observed)

# now, for each unique interaction, compute:
# - The proportion of islands where it was observed.
# - The average predicted probability (averaged over islands).
or_df_summary <- or_df_island %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )


or_df_filtered <- or_df %>%
  filter(itr == 1, original_links != 0) # filter existing interactions


# calculate proportion of islands in which each interaction occurs, rather than island pairs
or_df_filtered_self <- or_df_filtered %>% filter(train_layer == test_layer)

or_result_count <- or_df_filtered_self %>%
  # group by the interaction
  group_by(node_to, node_from) %>%
  summarise(
    # count unique test layers for this interaction
    n_test_layers = n_distinct(test_layer),
    .groups = "drop"
  ) %>%
  # calculate the proportion
  mutate(
    prop_test_layers = n_test_layers / n_distinct(df$test_layer)
  )


or_df_summary <- or_df_summary %>%
  left_join(or_result_count %>% select(node_to, node_from, prop_test_layers),
            by = c("node_to", "node_from")) %>%
  # replace avg_prop with the calculated proportion
  mutate(avg_prop_isl = if_else(is.na(prop_test_layers), 0, prop_test_layers)) %>%
  select(-prop_test_layers)  # remove helper column if not needed

# reset the levels for the species factors
or_df_summary$node_from <- factor(or_df_summary$node_from, levels = plant_order)
or_df_summary$node_to   <- factor(or_df_summary$node_to, levels = poll_order)

or_map_missing_links <- ggplot(or_df_summary, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop_isl)) +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = or_df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > or_best_discrete_threshold),
    aes(fill = avg_sigm_predicted),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2", 
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(or_map_missing_links)

pdf(
  file   = "results/paper_figs/map_missing_links_f1_as_cutoff.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(or_map_missing_links)
dev.off()     # close the file



