# ---- Subset analysis: comparing analyses quality ----
# here we compare prediction quality when using the entire species pools of A and P for prediction vs. prediction done based on subsetting only the species they share
# for reproducing Fig. 2a,b
## ---- load libraries ----
library(tidyverse)
library(emln)
library(grid)
library(ggpubr)
library(softImpute)
library(rstatix)
library(patchwork)

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

# set new layer names using the old ones
aggregated_df <- aggregated_df %>% 
  separate_wider_delim(layer_from, delim = "_", names = c("t", "l1", "l2"), cols_remove = FALSE) %>%
  mutate(island_id = paste0("layer_", as.numeric(l2)/2))  %>%
  mutate(layer_from = island_id, layer_to = island_id)%>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# View updated aggregated_df
print(aggregated_df)

# Total number of layers
num_layers <- length(unique(aggregated_df$layer_from))

results_file <- "results/predictions_island_scale.rds"

# read the prediction data if you already have it, and if not generate predictions

# Initialize a data frame to store combined results for all layer combinations
df_all <- data.frame()

if (file.exists(results_file)) {
  print("Existing results file found — reading the file and proceeding to analysis")
  
  df_all <- readRDS(results_file)
  print("finished loading prediction results")
  
} else { # or alternatively run the prediction pipeline
  # Loop through all combinations of layers_to_train and layer_to_predict
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
      
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
        
        sum(is.na(C))
        
        ### ---- b. + d. prediction with SVD and apply for all network combinations ----
        k_values <- c(2, 5, 10)
        lam0 <- lambda0(C)
        lambda_values <- c(1, 5, 50, 100, lam0)
        
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
            r <- implement_impute(C, k, lambda)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
            results <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }
        
        ### ---- save results for current k/lambda combination ----
        # After finishing the k/lambda loops, append the 'results' to 'df_all'
        # ---- (D) Append to df_all
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
  # df_all includes predictions for all combinations of islands, 50 iterations of links withholding and prediction for each combination
  # save the results
  saveRDS(df_all, file = results_file)
} 

# after reading or producing the results, filter these (important!):
df_all <- df_all %>% 
  filter(k == 2) %>% 
  filter(!(input_lambda %in%  c(1, 5, 50, 100)))

## ---- only shared species ----
set.seed(42)
# Initialize a data frame to store combined results for all layer combinations
df_shared_species <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
    
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
      
      sum(is.na(C))
      
      ### ---- b. + d. prediction with SVD ----
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
      # After finishing the k/lambda loops, append the 'results' to 'df_shared_species'
      # ---- (D) Append to df_shared_species
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
  
  list(
    best_threshold = best_threshold,
    result_summary = result_summary,
    metric_plot    = metric_plot
  )
}


res_all <- analyze_predictions(df_all, "All")
res_shared_species <- analyze_predictions(df_shared_species, "Shared species")

### ---- Fig. 2a,b: plot difference between analyses ----
df_all_combined <- bind_rows(
  res_all$result_summary             %>% mutate(dataset = "all"),
  res_shared_species$result_summary %>% mutate(dataset = "shared_species")
)
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
  shapiro_test(f1_score) # not normally distributed for all, use wilcoxon

# variance check
df_plot %>% levene_test(f1_score ~ dataset)

# nnse
# first normality check
df_plot %>%
  group_by(dataset) %>%
  shapiro_test(nnse) # nnse is normally distributed, use t-test

# variance check
df_plot %>% levene_test(nnse ~ dataset)

# make individual plots 
p_nnse <- ggplot(df_plot, aes(x = dataset, y = nnse, fill = dataset)) +
  geom_boxplot(notch = TRUE, alpha = 0.6) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "t.test",
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

# pdf(
#   file   = "p_nnse.pdf",
#   width  = 5,    # inches
#   height = 5,
#   family = "Helvetica"   # or another installed font
# )
# print(p_nnse)
# dev.off()     # close the file
# 

p_f1 <- ggplot(df_plot, aes(x = dataset, y = f1_score, fill = dataset)) +
  geom_boxplot(notch = TRUE, alpha = 0.6) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
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

# pdf(
#   file   = "p_f1.pdf",
#   width  = 5,    # inches
#   height = 5,
#   family = "Helvetica"   # or another installed font
# )
# print(p_f1)
# dev.off()     # close the file


# combine with a shared x-axis label
combined <- (p_nnse | p_f1) +       # side by side
  plot_layout(ncol = 2) &           # ensure two columns
  labs(x = "Network subset")     # shared x‐axis label

print(combined)

pdf(
  file   = "subset_analysis.pdf",
  width  = 7,    # inches
  height = 5,
  family = "Helvetica"   # or another installed font
)
print(combined)
dev.off()     # close the file

