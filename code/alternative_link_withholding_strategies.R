# degree based link withholding -
# for each node, save its degree. then there will be two degree values for each link - multiply them.
# then sample from the list of links, when the priority is in proportion to their degree:
# 1. positive proportionality - the higher the degree, the more likely to be sampled
# 2. negative proportionality - the higher the degree, the less likely to be sampled
# - the 3rd option is to sample uniformly at random, which is the same as not using degree information at all - the ones in the main script.
# output: figures 3a, 2c from the paper, for each option.


# load packages
source("code/common.R")
library(cowplot)
library(pROC)
library(PRROC)
library(ggnewscale)

# params ----
set.seed(42) # the answer to everything

# functions ----

# Find the threshold (on sigmoid-transformed predictions) that maximises the
# mean F0.5 score across all withheld observations in a results data frame.
find_optimal_threshold <- function(df, thresholds = seq(0, 1, by = 0.1)) {
  df_prepped <- df %>%
    filter(removed == 1) %>%
    mutate(
      predicted_prob  = sigmoid(predicted_values),
      original_binary = if_else(original_links > 0, 1, 0)
    )
  
  df_thresh <- df_prepped %>%
    tidyr::expand_grid(threshold = thresholds) %>%
    mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
    group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
    summarise(
      TP        = sum(original_binary == 1 & predicted_bin == 1),
      FP        = sum(original_binary == 0 & predicted_bin == 1),
      FN        = sum(original_binary == 1 & predicted_bin == 0),
      precision = TP / (TP + FP),
      recall    = TP / (TP + FN),
      f05_score = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
      .groups   = "drop"
    ) %>%
    group_by(emln_id, train_layer, test_layer, threshold) %>%
    summarise(f05_score = mean(f05_score, na.rm = TRUE), .groups = "drop") %>%
    group_by(threshold) %>%
    summarise(f05_score = mean(f05_score, na.rm = TRUE), .groups = "drop")
  
  best_row      <- df_thresh[which.max(df_thresh$f05_score), ]
  best_threshold <- best_row$threshold
  best_f05       <- best_row$f05_score
  
  message(sprintf("Optimal threshold: %.1f  (mean F0.5 = %.4f)", best_threshold, best_f05))
  return(best_threshold)
}

predict_with_degree_dependant_link_holdout <- function(aggregated_df, negative_degree_effect) {
  n_layers <- length(unique(aggregated_df$layer_from))
  
  # run predictions
  combined_results_degree_based <- NULL
  for (layers_to_train in 1:n_layers) {
    for (layer_to_predict in 1:n_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
      
      # calculate degree for each node in the training layers
      degree_from <- rowSums(P>0) # degree of node_from species
      degree_to <- colSums(P>0) # degree of node_to species
      
      ### ---- a. withhold links in P ----
      # map out the 0s and 1s in P
      ones_in_P <- which(P > 0, arr.ind = TRUE)
      zeros_in_P <- which(P == 0, arr.ind = TRUE)
      
      num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
      num_0_to_remove <- num_1_to_remove
      prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = T)
      
      
      # debug print
      print(paste("1 remove:", num_1_to_remove))
      print(paste("all 1   :", nrow(ones_in_P)))
      print(paste("0s to remove:", num_0_to_remove))
      print(paste("all zeros   :", nrow(zeros_in_P)))
      print(paste("prop of zeros removed   : ", prop_0_removed))
      
      
      # calculate the degree-based sampling probabilities for 1s and 0s - 
      # we use the product of the degrees of the two nodes as the sampling probability
      # For 1s
      degree_product <- degree_from[ones_in_P[, "row"]] * degree_to[ones_in_P[, "col"]]
      sampling_prob_1s <- degree_product / sum(degree_product)
      
      # For 0s
      degree_product_zeros <- degree_from[zeros_in_P[, "row"]] * degree_to[zeros_in_P[, "col"]]
      sampling_prob_0s <- degree_product_zeros / sum(degree_product_zeros)
      
      # if we want negative proportionality, we can use the inverse of the degree product
      if (negative_degree_effect) {
        sampling_prob_1s <- 1 / degree_product
        sampling_prob_1s <- sampling_prob_1s / sum(sampling_prob_1s) # normalize to sum to 1
        
        # for 0s
        sampling_prob_0s <- 1 / degree_product_zeros
        sampling_prob_0s <- sampling_prob_0s / sum(sampling_prob_0s)
      }
      
      # Randomly select zeros to withhold - bootstrapping
      bootstrapping_results <- NULL
      P_original <- P # save it for later
      
      # run the link withholding procedure with degree-based sampling
      for (i in 1:50) {
        # remove 1s
        remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove, prob = sampling_prob_1s), ]
        P[remove_indices] <- NA  # Set removed links to NA
        
        # sample 0s
        zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove, prob = sampling_prob_0s), ]
        P[zeros_to_remove_indices] <- NA
        
        ### ---- creating a combined matrix C ----
        # Combine A and P into a single matrix C with NAs representing missing data
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
        
        
        # Apply biScale to center matrices
        C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
        # save centers before overwriting (for back-transforming later)
        backtrans_vals <- c()
        backtrans_vals$row_centers <- attr(C, "biScale:row")$center      # named vector, length = nrow(C)
        backtrans_vals$col_centers <- attr(C, "biScale:column")$center   # named vector, length = ncol(C)
        sum(is.na(C))
        
        ### ---- b. + d. prediction with SVD and apply for all network combinations ----
        k_values <- c(2)
        lam0 <- lambda0(C)
        lambda_values <- c(lam0)
        
        # Initialize variables to store the best results
        results <- data.frame(k = integer(),
                              lambda = numeric(),
                              original_links = numeric(),
                              predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL
        
        # Loop over all combinations of k and lambda
        for (k in k_values) {
          for (lambda in lambda_values) {
            # imputation
            r <- implement_impute(C, k, lambda, P, remove_indices, 
                                  zeros_to_remove_indices, P_original, back_trans_values = backtrans_vals)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
            results <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }
        
        ### ---- save results for current k/lambda combination ----
        # After finishing the k/lambda loops, append the 'results' to 'combined_results'
        # ---- (D) Append to combined_results
        complete_edges_all <- rbind(results, not_removed_all)
        complete_edges_all$itr <- i
        bootstrapping_results <- rbind(bootstrapping_results, complete_edges_all)
        
        # reset P
        P <- P_original
      }
      
      combined_results_degree_based <- rbind(
        combined_results_degree_based,
        cbind(
          data.frame(
            emln_id = emln_id,
            train_layer = layers_to_train,
            test_layer = layer_to_predict,
            prop_ones_removed = prop_ones_to_remove,
            amount_of_removed_1 = num_1_to_remove,
            amount_of_removed_0 = num_0_to_remove,
            prop_0_removed = prop_0_removed
          ),
          bootstrapping_results
        )
      )
    }
  }
  
  return(combined_results_degree_based)
}

# run ----------
# load and mold data
aggregated_df <- load_and_mold_data_for_prediction(emln_id)

combined_results_degree_based_pos <- 
  predict_with_degree_dependant_link_holdout(aggregated_df, 
                                             negative_degree_effect = FALSE)
print("finished predicting with positive degree effect")

combined_results_degree_based_neg <-
  predict_with_degree_dependant_link_holdout(aggregated_df, 
                                             negative_degree_effect = TRUE)
print("finished predicting with negative degree effect")


# plot for results ----
# The plots: similar to figures 3a and 2c in the paper. 
# These will go to the supplementary.

## prepare data for plotting ----
# convert negatives to zeros
df_pos <- combined_results_degree_based_pos %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))
df_neg <- combined_results_degree_based_neg %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

prepare_results_to_plot <- function(df, threshold) {
  new_layer_names <- get_island_names_with_layer_indexes()
  
  result_summary <- df %>%
    filter(removed == 1) %>% 
    mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
    mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > threshold, 1, 0)) %>%
    mutate(original_binary = if_else(original_links > 0, 1, 0)) %>%
    group_by(emln_id, train_layer, test_layer, itr) %>% # start summarizing data
    summarise(
      TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
      FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
      TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
      FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
      specificity = TN / (TN + FP),
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      f05_score = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
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
      f05_score = mean(f05_score, na.rm = TRUE),
      balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
      nse  = mean(nse,  na.rm = TRUE),
      nnse = mean(nnse, na.rm = TRUE)
    ) %>%
    ungroup() %>%# Add island names to main table
    left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
    rename(train_layer_name = name) %>%
    left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
    rename(test_layer_name = name)
  
  return(result_summary)
}

threshold_pos <- find_optimal_threshold(df_pos)
threshold_neg <- find_optimal_threshold(df_neg)

result_summary_pos <- prepare_results_to_plot(df_pos, threshold_pos)
result_summary_neg <- prepare_results_to_plot(df_neg, threshold_neg)

## ---- Fig. S???: based on fig.2c - heatmap ----
plot_island_heatmap_with_degree_holdout <- function(result_summary, plot_title = NULL) {
  
  # Compute limits and midpoint dynamically
  lims <- range(result_summary$f05_score, na.rm = TRUE)
  mid_val <- mean(lims)
  
  island_heatmap_degree <- ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f05_score)) +
    
    # Base heatmap
    geom_tile(color = "black", linewidth = 0.1) +
    
    # Highlight diagonal
    geom_tile(
      data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
      color = "black", linewidth = 1.2
    ) +
    
    # Dynamic color scale
    scale_fill_gradient2(
      low = "lightsteelblue2",
      mid = "white",
      high = "salmon2",
      midpoint = mid_val,
      limits = lims,
      na.value = "gray"
    ) +
    
    labs(
      title = plot_title,
      x = "Added location",
      y = "Predicted location",
      fill = "F0.5 score"
    ) +
    
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    ) +
    
    coord_fixed() + tme
  
  return(island_heatmap_degree)
}

heatmap_pos <- 
  plot_island_heatmap_with_degree_holdout(result_summary_pos, 
                                          plot_title = "Positive degree effect")
heatmap_neg <-
  plot_island_heatmap_with_degree_holdout(result_summary_neg, 
                                          plot_title = "Negative degree effect")
print(heatmap_pos)
print(heatmap_neg)

fig_degree_holdout <- plot_grid(
  heatmap_pos,
  heatmap_neg,
  labels = c("(a)", "(b)"),
  ncol = 2,
  label_size = 17,
  label_x = 0.02,   # horizontal position (default ≈ 0)
  label_y = 0.78    # move labels closer to the plot
)

pdf(
  file   = "results/paper_figs/island_heatmap_degree_holdout_combined.pdf",
  width  = 11,    # inches
  height = 11,
  family = "Helvetica"   # or another installed font
)
print(fig_degree_holdout)
dev.off()     # close the file


## ---- Fig. S???: based on fig.3a ----
plot_link_prediction_map <- function(df, best_discrete_threshold, map_title) {
  # for each island and interaction, determine if the interaction was observed.
  # we use `any(original_links == 1)` so that if the interaction is observed in at least one iteration, we count it.
  df_island <- df %>%
    mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
    mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%
    group_by(node_from, node_to, island_id) %>%
    summarise(
      observed = as.integer(any(original_links != 0)),
      # For predicted_prob_sigm, you might take the average across iterations per island.
      island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
      .groups = "drop"
    )
  
  # now, for each unique interaction, compute:
  # - The proportion of islands where it was observed.
  # - The average predicted probability (averaged over islands).
  df_summary <- df_island %>%
    group_by(node_from, node_to) %>%
    summarise(
      avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
      avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
      n_islands = n(),  # number of islands contributing
      .groups = "drop"
    )
  
  # get degrees:
  df_filtered <- df %>%
    filter(itr == 1, original_links != 0) # filter existing interactions
  overall_poll_degree <- df_filtered %>% 
    group_by(node_to) %>% 
    summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")
  overall_plant_degree <- df_filtered %>% 
    group_by(node_from) %>% 
    summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")
  
  # now to plot:
  # here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of cases in which these interactions were observed.
  
  # order species by their degree
  df_summary <- df_summary %>% left_join(overall_poll_degree, by="node_to")
  df_summary <- df_summary %>% left_join(overall_plant_degree, by="node_from")
  
  # determine the order of species in the plot based on their degree
  plant_order <- df_summary %>%
    distinct(node_from, overall_plant_degree) %>%
    arrange(desc(overall_plant_degree)) %>%
    pull(node_from)
  
  poll_order <- df_summary %>%
    distinct(node_to, overall_poll_degree) %>%
    arrange(desc(overall_poll_degree)) %>%
    pull(node_to)
  
  # calculate proportion of islands in which each interaction occurs, rather than island pairs
  df_filtered_self <- df_filtered %>% filter(train_layer == test_layer)
  
  result_count <- df_filtered_self %>%
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
  
  df_summary <- df_summary %>%
    left_join(result_count %>% select(node_to, node_from, prop_test_layers),
              by = c("node_to", "node_from")) %>%
    # replace avg_prop with the calculated proportion
    mutate(avg_prop_isl = if_else(is.na(prop_test_layers), 0, prop_test_layers)) %>%
    select(-prop_test_layers)  # remove helper column if not needed
  
  # reset the levels for the species factors
  df_summary$node_from <- factor(df_summary$node_from, levels = plant_order)
  df_summary$node_to   <- factor(df_summary$node_to, levels = poll_order)
  
  map_missing_links_degree_based <- ggplot(df_summary, aes(x = node_to, y = node_from)) +
    # First layer: background heatmap for proportion observed (blue gradient)
    geom_tile(aes(fill = avg_prop_isl)) +
    scale_fill_gradient(low = "white", high = "steelblue", 
                        name = "Observed links:\nproportion\nof islands\nobserved",
                        breaks = seq(0, 1, 0.2)) +
    
    # Reset fill scale so the next layer can have its own gradient
    new_scale_fill() +
    
    # Second layer: overlay only cells that were never observed but have high predicted value
    geom_tile(
      data = df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
      aes(fill = avg_sigm_predicted),
      alpha = 0.6
    ) +
    scale_fill_gradient(low = "tan1", high = "tomato2", 
                        name = "Predicted links:\naverage predicted\nprobability",
                        breaks = seq(0, 1, 0.1)) +
    
    # Final adjustments
    theme_minimal() +
    labs(title = map_title, x = "Pollinator", y = "Plant") +
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
  
  return(map_missing_links_degree_based)
}

map_links_degree_pos <-
  plot_link_prediction_map(df_pos, threshold_pos,
                           map_title = "Degree-based withholding: positive degree effect")
map_links_degree_neg <-
  plot_link_prediction_map(df_neg, threshold_neg,
                           map_title = "Degree-based withholding: negative degree effect")

print(map_links_degree_pos)
print(map_links_degree_neg)

pdf(
  file   = "results/paper_figs/map_missing_links_degree_based_holdout.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(map_links_degree_pos)
print(map_links_degree_neg)
dev.off()     # close the file


# ---- class-imbalance dependent withholding ----
predict_with_class_imbalance_link_holdout <- function(aggregated_df) {
  n_layers <- length(unique(aggregated_df$layer_from))
  
  # run predictions
  combined_results_class_imbalance <- NULL
  for (layers_to_train in 1:n_layers) {
    for (layer_to_predict in 1:n_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
      
      ### ---- a. withhold links in P ----
      # map out the 0s and 1s in P
      ones_in_P <- which(P > 0, arr.ind = TRUE)
      zeros_in_P <- which(P == 0, arr.ind = TRUE)
      
      # set proportions to remove
      prop_ones_to_remove  <- 0.2
      prop_zeros_to_remove <- 0.2
      
      # sample links and non-links
      num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
      num_0_to_remove <- floor(sum(P == 0, na.rm = TRUE) * prop_zeros_to_remove)  # Number of non-links to remove
      prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = TRUE)
      
      # debug print
      print(paste("1 remove:", num_1_to_remove))
      print(paste("all 1   :", nrow(ones_in_P)))
      print(paste("0s to remove:", num_0_to_remove))
      print(paste("all zeros   :", nrow(zeros_in_P)))
      print(paste("prop of zeros removed   : ", prop_0_removed))
      
      # Randomly select zeros to withhold - bootstrapping
      bootstrapping_results <- NULL
      P_original <- P # save it for later
      
      # run the link withholding procedure with degree-based sampling
      for (i in 1:50) {
        # remove 1s
        remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
        P[remove_indices] <- NA  # Set removed links to NA
        
        # sample 0s
        zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
        P[zeros_to_remove_indices] <- NA
        
        ### ---- creating a combined matrix C ----
        # Combine A and P into a single matrix C with NAs representing missing data
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
        
        
        # Apply biScale to center matrices
        C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
        # save centers before overwriting (for back-transforming later)
        backtrans_vals <- c()
        backtrans_vals$row_centers <- attr(C, "biScale:row")$center      # named vector, length = nrow(C)
        backtrans_vals$col_centers <- attr(C, "biScale:column")$center   # named vector, length = ncol(C)
        sum(is.na(C))
        
        ### ---- b. + d. prediction with SVD and apply for all network combinations ----
        k_values <- c(2)
        lam0 <- lambda0(C)
        lambda_values <- c(lam0)
        
        # Initialize variables to store the best results
        results <- data.frame(k = integer(),
                              lambda = numeric(),
                              original_links = numeric(),
                              predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL
        
        # Loop over all combinations of k and lambda
        for (k in k_values) {
          for (lambda in lambda_values) {
            # imputation
            r <- implement_impute(C, k, lambda, P, remove_indices, 
                                  zeros_to_remove_indices, P_original, back_trans_values = backtrans_vals)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
            results <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }
        
        ### ---- save results for current k/lambda combination ----
        # After finishing the k/lambda loops, append the 'results' to 'combined_results'
        # ---- (D) Append to combined_results
        complete_edges_all <- rbind(results, not_removed_all)
        complete_edges_all$itr <- i
        bootstrapping_results <- rbind(bootstrapping_results, complete_edges_all)
        
        # reset P
        P <- P_original
      }
      
      combined_results_class_imbalance <- rbind(
        combined_results_class_imbalance,
        cbind(
          data.frame(
            emln_id = emln_id,
            train_layer = layers_to_train,
            test_layer = layer_to_predict,
            prop_ones_removed = prop_ones_to_remove,
            amount_of_removed_1 = num_1_to_remove,
            amount_of_removed_0 = num_0_to_remove,
            prop_0_removed = prop_0_removed
          ),
          bootstrapping_results
        )
      )
    }
  }
  
  return(combined_results_class_imbalance)
}


combined_results_class_imbalance <-
  predict_with_class_imbalance_link_holdout(aggregated_df)
print("finished predicting with class imbalance")

df_im <- combined_results_class_imbalance %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))  # convert the predicted values to probability values in the interval (0, 1) using the logistic function


threshold_im <- find_optimal_threshold(df_im)
result_summary_imbalance <- prepare_results_to_plot(df_im, threshold_im)

## ---- plot roc and pr curves for class imbalance ----
# positive class prevalence
prev_pos <- mean(df_im$original_binary == 1, na.rm = TRUE)

# ROC object
roc_obj <- roc(
  response = df_im$original_binary,
  predictor = df_im$predicted_prob_sigm,
  quiet = TRUE,
  na.rm = TRUE,
  levels = c(0, 1),
  direction = "<"
)

auc_roc_value <- as.numeric(auc(roc_obj))

# point for chosen threshold
threshold <- threshold_im

coords_df <- coords(
  roc_obj,
  x = threshold,
  input = "threshold",
  ret = c("specificity", "sensitivity")
)


# full ROC dataframe
roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

roc_point <- data.frame(
  FPR = 1 - coords_df[["specificity"]],
  TPR = coords_df[["sensitivity"]]
)

roc_plot <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth = 1.2, color = "lightsteelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "salmon") +
  #geom_point(data = roc_point, aes(FPR, TPR), size = 3, color = "black") +
  labs(
    title = "ROC curve",
    subtitle = paste0("AUC = ", round(auc_roc_value, 3)),
    x = "False positive rate (1 - specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_classic(base_size = 14) +
  tme

# pr
df_im <- df_im %>% filter(removed == 1)

pos_scores <- df_im$predicted_prob_sigm[df_im$original_binary == 1]
neg_scores <- df_im$predicted_prob_sigm[df_im$original_binary == 0]

pr_obj <- pr.curve(
  scores.class0 = pos_scores,
  scores.class1 = neg_scores,
  curve = TRUE
)

auc_pr_value <- pr_obj$auc.integral

pr_df <- data.frame(
  Recall = pr_obj$curve[, 1],
  Precision = pr_obj$curve[, 2]
)

pr_plot <- ggplot(pr_df, aes(Recall, Precision)) +
  geom_line(linewidth = 1.2, color = "lightsteelblue") +
  geom_hline(yintercept = prev_pos, linetype = "dashed", color = "salmon") +
  labs(
    title = "Precision-Recall curve",
    subtitle = paste0(
      "PR-AUC = ", round(auc_pr_value, 3),
      " | baseline = ", round(prev_pos, 3)
    ),
    x = "Recall",
    y = "Precision"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_classic(base_size = 14) +
  tme

fig_imbalance_pr_roc <- plot_grid(
  roc_plot,
  pr_plot,
  labels = c("(a)", "(b)"),
  ncol = 2,
  label_size = 17,
  label_x = 0.02,   # horizontal position (default ≈ 0)
  label_y = 1.1    # move labels closer to the plot
)

fig_imbalance_pr_roc

pdf(
  file   = "results/paper_figs/fig_imbalance_pr_roc.pdf",
  width  = 11,    # inches
  height = 11,
  family = "Helvetica"   # or another installed font
)
print(fig_imbalance_pr_roc)
dev.off()     # close the file


# roc heatmap

df_eval_im <- df_im %>%
  filter(removed == 1) %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    # ROC-AUC (coerce to numeric!)
    auc_roc = tryCatch({
      roc_obj <- roc(response = original_binary,
                     predictor = predicted_prob_sigm,
                     quiet = TRUE, na.rm = TRUE,
                     levels = c(0,1), direction = "<")
      as.numeric(auc(roc_obj))   # <-- important
    }, error = function(e) NA_real_),
    
    # PR-AUC (guard against all-one-class cases)
    auc_pr = tryCatch({
      pos <- predicted_prob_sigm[original_binary == 1]
      neg <- predicted_prob_sigm[original_binary == 0]
      if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
      pr_obj <- pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)
      pr_obj$auc.integral
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

df_eval_summary_im <- df_eval_im %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    auc_roc_mean = mean(auc_roc, na.rm = TRUE),
    auc_roc_sd   = sd(auc_roc,   na.rm = TRUE),
    auc_pr_mean  = mean(auc_pr,  na.rm = TRUE),
    auc_pr_sd    = sd(auc_pr,    na.rm = TRUE),
    .groups = "drop"
  )

# add layer names
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# create a new grouped tibble for island names
new_layer_names <- net_name %>%
  mutate(group_id = (layer_id + 1) %/% 2) %>%  # Group pairs into 1, 2, 3...
  group_by(group_id) %>%
  summarise(name = gsub(" site.*", "", first(name)), .groups = "drop")  # Keep only location name

# Add names to main table
df_eval_summary_im <- df_eval_summary_im %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)


# roc
imbalance_heatmap_roc <- 
  ggplot(df_eval_summary_im, aes(x = train_layer_name, y = test_layer_name, fill = auc_roc_mean)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = df_eval_summary_im[df_eval_summary_im$train_layer == df_eval_summary_im$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "thistle", 
                       midpoint = 0.69, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added location", y = "Predicted location", fill = "ROC-AUC") +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(imbalance_heatmap_roc)

