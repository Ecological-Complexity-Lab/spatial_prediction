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

# ------------- parsing arguments -----------
# read args given in command line:
if (length(commandArgs(trailingOnly=TRUE))==0) { # make sure we have commands
  stop('No arguments were found!') # the script will not run without arguments
} else {
  args <- commandArgs(trailingOnly=TRUE)
  emln_id <- as.numeric(args[1])
  is_binary <- as.numeric(args[2])
  
}

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

implement_impute <- function(C, k, lambda) {
  # Apply softImpute
  
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  
  # Debias the fit to remove regularization effects
  # fit <- deBias(C, fit)
  
  # Reconstruct the matrix
  C_reconstructed <- softImpute::complete(C, fit, unscale = TRUE)
  
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

## ---- parameters ----
emln_id <- 60
layers_to_train <- 1
layer_to_predict <- 7
prop_ones_to_remove <- 0.2
# prop_zeros_to_remove <- 0.2
n_sim <- 50
is_binary <- 0
set.seed(42)

## ---- create a folder for the results ----
# setwd("~/softimpute/results")
# 
# if (!dir.exists("results_net_60_weighted_50_itr")) {
#   dir.create("results_net_60_weighted_50_itr", recursive = TRUE)
# }

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

# Save the mapping to CSV
# write_csv(layer_mapping, "layer_mapping.csv")

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
combined_results <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    # # # crop to make toy example matrices
    #A <- A[1:8, 1:7]
    #P <- P[1:8, 1:7]
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    if (is_binary == 1) {
      A[A > 0] <- 1  # Make binary
      P[P > 0] <- 1  # Make binary
    }
    
    ### ---- remove some links in P ----
    # map out the 0s and 1s in P
    num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
    ones_in_P <- which(P > 0, arr.ind = TRUE)
    
    num_0_to_remove <- num_1_to_remove
    prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = T)
    #num_0_to_remove <- floor((nrow(P)*ncol(P)-sum(P>0, na.rm = T))*prop_zeros_to_remove)  # Number of links to remove
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    
    # debug print
    print(paste("1 remove:", num_1_to_remove))
    print(paste("all 1   :", nrow(ones_in_P)))
    print(paste("0s to remove:", num_0_to_remove))
    print(paste("all zeros   :", nrow(zeros_in_P)))
    print(paste("prop of zeros removed   : ", prop_0_removed))
    
    # # remove 1s
    # remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
    # P[remove_indices] <- NA  # Set removed links to NA
    # P_no_1 <- P # save it for later
    
    bootstrapping_results <- NULL
    P_original <- P # save it for later
    
    # Randomly select zeros to remove - bootstrapping
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
      
      if (is_binary == 1) {
        C[C>0] <- 1 # Make binary
      }
      
      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
      
      sum(is.na(C))
      # might need to convert C into a binary matrix
      
      ### ---- transfer learning with SVD ----
      # Define the grid of k and lambda values to search over
      #k_values <- c(2, 3, 4, 5, 10, 15, 20)            # Adjust as needed
      # lambda_values <- c(0, 0.001, 0.01, 0.05, 0.1)  # Adjust as needed
      
      k_values <- c(2)
      
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
    
    combined_results <- rbind(
      combined_results,
      cbind(
        data.frame(
          emln_id = emln_id,
          train_layer = layers_to_train,
          test_layer = layer_to_predict,
          prop_ones_removed = prop_ones_to_remove,
          amount_of_removed_1 = num_1_to_remove,
          #prop_zeros_removed = prop_zeros_to_remove,
          amount_of_removed_0 = num_0_to_remove,
          prop_0_removed = prop_0_removed
        ),
        bootstrapping_results
      )
    )
  }
}

# View combined results
summary(combined_results)

# Save the combined results dataframe to a CSV file
#output_name <- paste0("binary_equal_0_1_removal_scaling_island_",emln_id,"_",is_binary,".csv")
output_name <- paste0("_unscale_weighted__scaled_island_net_",emln_id,"_",n_sim,"_itr.csv")
write.csv(combined_results, file = output_name, row.names = FALSE)
#write.csv(df, file = "duplicate_check.csv", row.names = FALSE)
