# ---- Subset analysis: comparing analyses quality ----
# here we compare prediction quality when using the entire species pools of A and P for prediction vs. prediction done based on subsetting only the species they share
## ---- load libraries ----
library(tidyverse)
library(pROC)
library(PRROC)
library(emln)
library(grid)
library(ggpubr)

## ---- parameters ----
emln_id <- 60 # Canary Islands pollination system from Trøjelsgaard et al. 2015
prop_ones_to_remove <- 0.2 # proportion of existing links to withhold
n_sim <- 50 # number of random link withholding and prediction iterations
set.seed(42) # the answer to everything

## ---- functions ----
# building matrices for combining matrices, calculating network size and density
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

# predict links using softImpute
implement_impute <- function(C, k, lambda) {
  # Apply softImpute
  
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  
  # Debias the fit to remove regularization effects
  # fit <- deBias(C, fit)
  
  # Reconstruct the matrix
  C_reconstructed <- softImpute::complete(C, fit)
  
  # Extract the reconstructed P matrix from C_reconstructed
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
  # Combine indices of removed ones and zeros
  if (is.null(dim(remove_indices))) { # handles when remove_indices has only one row
    test_indices <- rbind(
      data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  } else {
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = rep(1, nrow(remove_indices))),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  }
  
  # Get the row and column names of the test links
  test_rows <- rownames(P)[test_indices$row]  # These are the "node_to"
  test_cols <- colnames(P)[test_indices$col]  # These are the "node_from"
  
  # Actual labels and predictions
  original_links <- P_original[cbind(test_rows, test_cols)]
  predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
  
  # Store the results with node information
  results <- data.frame(k = k,
                        lambda = lambda,
                        original_links = original_links,
                        predicted_values = predicted_values,
                        node_to = test_rows,
                        node_from = test_cols,
                        removed = 1 # mark these links as removed
  )
  
  # ---- (A) Enumerate ALL edges in P
  all_edges <- expand.grid(
    node_to = rownames(P),
    node_from = colnames(P),
    k = k,
    lambda = lambda,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Fill in the original link values from P_original
  all_edges$original_links <- mapply(
    function(r, c) P_original[r, c],
    all_edges$node_to,
    all_edges$node_from
  )
  
  # Helper data frame of removed edges
  removed_edges_idx <- data.frame(
    node_to = test_rows,
    node_from = test_cols,
    stringsAsFactors = FALSE
  )
  
  # ---- (B) Subset edges NOT removed
  not_removed <- all_edges[
    !paste(all_edges$node_to, all_edges$node_from) %in%
      paste(removed_edges_idx$node_to, removed_edges_idx$node_from), # all node pairs that are not in th removed list
  ]
  not_removed$removed <- 0
  not_removed$predicted_values <- NA
  
  # ---- (C) Combine removed + not removed
  not_removed$k <- k
  not_removed$lambda <- lambda
  
  list_results <- list(results=results, not_removed=not_removed)
  
  return(list_results)
}

# transformation

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Function to plot ROC curve with ggplot2
plot_roc_curve <- function(true_labels, predicted_scores) {
  # Create the ROC object and compute AUC
  roc_obj <- roc(response = true_labels, predictor = predicted_scores)
  auc_val <- auc(roc_obj)
  
  # Build a data frame from the ROC object for ggplot2
  df_roc <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  
  # Generate the ROC plot
  p <- ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # ROC curve line
    geom_abline(intercept = 0, slope = 1,                       # Diagonal line (random classifier)
                linetype = "dashed", color = "salmon") +
    labs(title = paste("ROC curve (AUC =", round(auc_val, 2), ")"),
         x = "False positive rate", y = "True positive rate") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

# Function to plot PR curve with ggplot2
plot_pr_curve <- function(true_labels, predicted_scores) {
  # Separate scores by class: positive (true label==1) and negative (true label==0)
  scores_pos <- predicted_scores[true_labels == 1]
  scores_neg <- predicted_scores[true_labels == 0]
  
  # Create the PR curve object; curve=TRUE returns the full curve data
  pr_obj <- pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = TRUE)
  
  # Calculate the positive class ratio for the random baseline line
  pos_ratio <- sum(true_labels == 1) / length(true_labels)
  
  # Convert the curve matrix to a data frame and set column names
  df_pr <- as.data.frame(pr_obj$curve)
  colnames(df_pr) <- c("recall", "precision", "threshold")
  
  # Generate the PR plot
  p <- ggplot(df_pr, aes(x = recall, y = precision)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # PR curve line
    geom_hline(yintercept = pos_ratio,                          # Baseline: random classifier performance
               linetype = "dashed", color = "salmon") +
    labs(title = paste("PR curve (AUC =", round(pr_obj$auc.integral, 2), ")"),
         x = "Recall", y = "Precision") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- load/produce data ----
### ---- all species ----
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

# Initialize a data frame to store combined results for all layer combinations
df_all <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    ### ---- a. withhold links in P ----
    # map out the 0s and 1s in P
    num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
    ones_in_P <- which(P > 0, arr.ind = TRUE)
    
    num_0_to_remove <- num_1_to_remove
    prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = T)
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    
    # debug print
    print(paste("1 remove:", num_1_to_remove))
    print(paste("all 1   :", nrow(ones_in_P)))
    print(paste("0s to remove:", num_0_to_remove))
    print(paste("all zeros   :", nrow(zeros_in_P)))
    print(paste("prop of zeros removed   : ", prop_0_removed))
    
    # Randomly select zeros to withhold - bootstrapping
    bootstrapping_results <- NULL
    P_original <- P # save it for later
    
    for (i in 1:n_sim) {
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
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
      
      
      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
      
      sum(is.na(C))
      
      ### ---- b. prediction with SVD ----
      k_values <- 2
      lam0 <- lambda0(C)
      lambda_values <- c(lam0)
      
      # Initialize variables to store the best results
      results <- data.frame(k = integer(),
                            lambda = numeric(),
                            original_links = numeric(),
                            predicted_values = numeric())
      not_removed_all <- NULL
      
      # Loop over all combinations of k and lambda
      for (k in k_values) {
        for (lambda in lambda_values) {
          # imputation
          r <- implement_impute(C, k, lambda)
          
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
    
    df_all <- rbind(
      df_all,
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

#df_all <- read.csv('weighted__scaled_island_net_60_100_itr.csv')

### ---- only shared species ----
set.seed(42)
# Initialize a data frame to store combined results for all layer combinations
df_shared_species <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    # 1) find the species they have in common
    shared_pollinators <- intersect(rownames(A), rownames(P))
    shared_plants      <- intersect(colnames(A), colnames(P))
    
    # 2) subset both matrices to exactly those shared species
    A <- A[shared_pollinators, shared_plants, drop = FALSE]
    P <- P[shared_pollinators, shared_plants, drop = FALSE]
    
    # now A and P have identical dimnames, and you can safely combine them:
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    ### ---- a. withhold links in P ----
    # map out the 0s and 1s in P
    num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
    ones_in_P <- which(P > 0, arr.ind = TRUE)
    
    num_0_to_remove <- num_1_to_remove
    prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = T)
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    
    # debug print
    print(paste("1 remove:", num_1_to_remove))
    print(paste("all 1   :", nrow(ones_in_P)))
    print(paste("0s to remove:", num_0_to_remove))
    print(paste("all zeros   :", nrow(zeros_in_P)))
    print(paste("prop of zeros removed   : ", prop_0_removed))
    
    if (num_1_to_remove < 2 | num_0_to_remove < 2) {
      print("1s or 0s too low; moving to the next layer pair")
      next
    } # because now our matrices are smaller
    
    # Randomly select zeros to withhold - bootstrapping
    bootstrapping_results <- NULL
    P_original <- P # save it for later
    
    for (i in 1:n_sim) {
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
      C <- matrix(0,
                  nrow = length(shared_pollinators),
                  ncol = length(shared_plants),
                  dimnames = list(shared_pollinators, shared_plants))
      
      # Place A into C
      C[rownames(A), colnames(A)] <- A
      
      # Place P into C
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
      
      
      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
      
      sum(is.na(C))
      
      ### ---- b. prediction with SVD ----
      k_values <- 2
      lam0 <- lambda0(C)
      lambda_values <- c(lam0)
      
      # Initialize variables to store the best results
      results <- data.frame(k = integer(),
                            lambda = numeric(),
                            original_links = numeric(),
                            predicted_values = numeric())
      not_removed_all <- NULL
      
      # Loop over all combinations of k and lambda
      for (k in k_values) {
        for (lambda in lambda_values) {
          # imputation
          r <- implement_impute(C, k, lambda)
          
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
    
    df_shared_species <- rbind(
      df_shared_species,
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

# or read a file:
#df_shared_species <- read.csv('weighted__scaled_island_net_60_100_itr_shared_species.csv')

# convert negative values to zeros
df_all <- df_all %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

df_shared_species <- df_shared_species %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

## ---- ok go ----
### ---- evaluation ----
analyze_predictions <- function(df, name = "Dataset") {
  thresholds <- seq(0, 1, by = 0.1)
  
  df_prepped <- df %>%
    filter(removed == 1) %>%
    mutate(
      predicted_prob   = sigmoid(predicted_values),
      original_binary  = if_else(original_links > 0, 1, 0)
    )
  
  df_thresh <- df_prepped %>%
    tidyr::expand_grid(threshold = thresholds) %>%
    mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
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
      rmse = sqrt(mse),
      .groups = "drop"
    ) %>%
    group_by(emln_id, train_layer, test_layer, threshold) %>%
    summarise(across(c(TP:rmse), mean, na.rm = TRUE), .groups = "drop")
  
  df_avg <- df_thresh %>%
    group_by(threshold) %>%
    summarise(across(
      c(specificity, precision, recall,
        f1_score, balanced_accuracy, mcc),
      mean, na.rm = TRUE
    )) %>%
    pivot_longer(-threshold,
                 names_to  = "metric",
                 values_to = "value")
  
  # Plot metric curves
  metric_plot <- ggplot(df_avg, aes(threshold, value, color = metric)) +
    geom_line(size = 1) +
    labs(title = paste(name, "- Average metrics"), x = "Threshold", y = "Value") +
    theme_minimal()
  
  df_wide <- df_avg %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(absdiff = abs(f1_score - balanced_accuracy)) %>%
    slice_min(absdiff, n = 1)
  
  best_threshold <- df_wide$threshold
  
  # Evaluate at optimal threshold
  df_removed <- df %>%
    filter(removed == 1) %>%
    mutate(
      predicted_prob_sigm = sigmoid(predicted_values),
      predicted_bin_sigm  = if_else(predicted_prob_sigm > best_threshold, 1, 0),
      original_binary      = if_else(original_links > 0, 1, 0)
    )
  
  result_summary <- df_removed %>%
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
      mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
      mse = mean((predicted_values - original_links)^2, na.rm = TRUE),
      rmse = sqrt(mse),
      nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
        sum((original_links   - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
      nnse = 1 / (2 - nse),
      .groups = "drop"
    ) %>%
    group_by(emln_id, train_layer, test_layer) %>%
    summarise(across(TP:nnse, mean, na.rm = TRUE), .groups = "drop")
  
  # ROC and PR Curves
  auc_val <- plot_roc_curve(df_removed$original_binary, df_removed$predicted_values)
  plot_pr_curve(df_removed$original_binary, df_removed$predicted_values)
  
  list(
    best_threshold = best_threshold,
    result_summary = result_summary,
    metric_plot    = metric_plot,
    auc            = auc_val
  )
}


res_all <- analyze_predictions(df_all, "All")
res_shared_species <- analyze_predictions(df_shared_species, "Shared species")

### ---- roc curves ----
# 1. Put your results into a named list
roc_data_list <- list(
  All                 = df_all,
  Shared_species      = df_shared_species
)

# 2. For each subset, compute a pROC::roc object and turn it into a data.frame
roc_df <- purrr::imap_dfr(roc_data_list, function(df, subset_name) {
  df2 <- df %>%
    filter(removed == 1) %>%
    mutate(
      truth = as.numeric(original_links > 0),
      prob  = sigmoid(predicted_values)
    )
  roc_obj <- roc(df2$truth, df2$prob)
  
  # pull out the sensitivities / specificities
  tibble(
    subset     = subset_name,
    fpr        = 1 - roc_obj$specificities,
    tpr        =    roc_obj$sensitivities,
    threshold  =    roc_obj$thresholds,
    auc        = as.numeric(roc_obj$auc)  # same for every row
  )
})

# 3. Extract one AUC‐per‐subset for annotation
auc_labels <- roc_df %>%
  select(subset, auc) %>%
  distinct() %>%
  mutate(label = paste0("AUC = ", round(auc, 3)),
         x = 0.6,  # x/y are positions in FPR‐TPR space
         y = 0.2)

# 4. Plot with facets
# 4. plot, with custom colors + labeller
#    change these to whatever you like:
line_color   <- "lightsteelblue"
random_color <- "salmon"

ggplot(roc_df, aes(x = fpr, y = tpr)) +
  # your ROC line
  geom_line(color = line_color, size = 1) +
  # diagonal random‐guess line
  geom_abline(
    intercept = 0, slope = 1,
    linetype  = "dashed",
    color     = random_color
  ) +
  # facet and relabel
  facet_wrap(
    ~ subset, ncol = 2,
    labeller = labeller(subset = c(
      All                = "All",
      Shared_Plants      = "Shared plants",
      Shared_Pollinators = "Shared pollinators",
      Shared_Species     = "Shared species"
    ))
  ) +
  # annotate each panel with its AUC
  geom_text(
    data        = auc_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 0
  ) +
  coord_equal() +                     # square 0–1 axes
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "False Positive Rate (1 − specificity)",
    y = "True Positive Rate (sensitivity)"
  ) +
  theme_minimal(base_size = 14) + tme

# pr curves

pr_data_list <- list(
  All                 = df_all,
  Shared_Species      = df_shared_species
)

pr_df <- imap_dfr(pr_data_list, function(df, subset_name) {
  df2 <- df %>%
    filter(removed == 1) %>%
    mutate(
      truth = as.numeric(original_links > 0),
      prob  = sigmoid(predicted_values)
    )
  
  pr_obj <- pr.curve(
    scores.class0 = df2$prob[df2$truth == 1],
    scores.class1 = df2$prob[df2$truth == 0],
    curve = TRUE
  )
  
  tibble(
    subset    = subset_name,
    recall    = pr_obj$curve[, 1],
    precision = pr_obj$curve[, 2],
    threshold = pr_obj$curve[, 3],
    aucpr     = pr_obj$auc.integral
  )
})

# 3. extract one AUPRC per subset for labeling
aucpr_labels <- pr_df %>%
  distinct(subset, aucpr) %>%
  mutate(
    label = paste0("AUPRC = ", round(aucpr, 3)),
    x = 0.2,  # adjust as needed
    y = 0.8
  )

# 4. choose your line color (reuse from ROC if you like)
line_color <- "lightsteelblue"

# 5. plot side-by-side PR curves
ggplot(pr_df, aes(x = recall, y = precision)) +
  # PR curve
  geom_line(color = line_color, size = 1) +
  # facet panels with pretty titles
  facet_wrap(
    ~ subset, ncol = 2,
    labeller = labeller(subset = c(
      All                = "All",
      Shared_Species     = "Shared species"
    ))
  ) +
  # annotate each panel with its AUPRC
  geom_text(
    data        = aucpr_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 0
  ) +
  labs(
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 14) + tme

### ---- plot difference between analyses ----
df_all_combined <- bind_rows(
  res_all$result_summary             %>% mutate(dataset = "all"),
  res_shared_species$result_summary %>% mutate(dataset = "shared_species")
)
# 
# # plot
# # overall difference
# # For nNSE
# ggplot(df_all_combined, aes(x = dataset, y = nnse)) +
#   geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   stat_compare_means(method = "anova", label = "p.format",
#                      label.y = max(df_all_combined$nnse, na.rm = TRUE) * 1.05) +
#   labs(title = "nNSE across datasets", x = "Dataset", y = "nNSE") +
#   theme_minimal() + tme
# 
# # For F1
# ggplot(df_all_combined, aes(x = dataset, y = f1_score)) +
#   geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
#   geom_jitter(width = 0.15, alpha = 0.5) +
#   stat_compare_means(method = "anova", label = "p.format",
#                      label.y = max(df_all_combined$f1_score, na.rm = TRUE) * 1.05) +
#   labs(x = "Dataset", y = "F1 Score") +
#   theme_minimal() + tme
# 
# comparisons <- list(
#   c("all","shared_species"), c("all","shared_plants"), c("all","shared_pollinators"),
#   c("shared_species","shared_plants"), c("shared_species","shared_pollinators"), c("shared_plants","shared_pollinators")
# )
# 
# # pairwise nnse
# # ggplot(df_all_combined, aes(x = dataset, y = nnse)) +
# #   geom_boxplot(notch = TRUE, alpha = 0.7) +
# #   stat_compare_means(
# #     comparisons = comparisons,
# #     method      = "wilcox.test",
# #     label       = "p.signif"
# #   ) +
# #   labs(x = "Network subset", y = "nNSE") +
# #   theme_minimal() + tme
# 
# 1) Tidy up + rename the groups
df_plot <- df_all_combined %>%
  mutate(
    dataset = factor(
      dataset,
      # original names in the data
      levels = c("all", "shared_species"),
      # how they should appear on the plot
      labels = c("All", "Shared species")
    )
  )

# 2) Define a custom palette
my_palette <- c(
  "All"                = "rosybrown2",
  "Shared species"     = "thistle"
)

# define comparison and palette:
comparisons <- list(c("All", "Shared species"))

# do we need welch/wilcoxon?
# first normality check
df_plot %>%
  group_by(dataset) %>%
  shapiro_test(f1_score)

# variance check
df_plot %>% levene_test(f1_score ~ dataset)

# nnse
# first normality check
df_plot %>%
  group_by(dataset) %>%
  shapiro_test(nnse)

# variance check
df_plot %>% levene_test(nnse ~ dataset)

# make individual plots 
p_nnse <- ggplot(df_plot, aes(x = dataset, y = nnse, fill = dataset)) +
  geom_boxplot(notch = TRUE, alpha = 0.6) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.format",
    tip.length  = 0.01
  ) +
  scale_fill_manual(values = my_palette, guide = FALSE) +
  labs(y = "NNSE") +
  theme_minimal() +
  tme +  # your custom theme
  theme(
    axis.text.x         = element_text(angle = 25, hjust = 1),
    panel.grid.major.x  = element_blank(),
    axis.title.x        = element_blank()  # we'll add a shared x‐label later
  )

p_f1 <- ggplot(df_plot, aes(x = dataset, y = f1_score, fill = dataset)) +
  geom_boxplot(notch = TRUE, alpha = 0.6) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "t.test",
    label       = "p.format",
    tip.length  = 0.01
  ) +
  scale_fill_manual(values = my_palette, guide = FALSE) +
  labs(y = "F1 score") +
  theme_minimal() +
  tme +
  theme(
    axis.text.x         = element_text(angle = 25, hjust = 1),
    panel.grid.major.x  = element_blank(),
    axis.title.x        = element_blank()
  )

# combine with a shared x-axis label
combined <- (p_nnse | p_f1) +       # side by side
  plot_layout(ncol = 2) &           # ensure two columns
  labs(x = "Network subset")     # shared x‐axis label

print(combined)

# pdf(
#   file   = "subset_analysis.pdf",
#   width  = 7,    # inches
#   height = 5,
#   family = "Helvetica"   # or another installed font
# )
# print(combined)
# dev.off()     # close the file

