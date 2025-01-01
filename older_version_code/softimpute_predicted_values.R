# ---- softImpute for predicting removed links in empirical networks - all layer combinations ----
# including removal of zeros and ones. no bootstrapping yet
# in this code we remove zeros and ones and get the predicted values (not labels but values based on the svd from the softImpute fit function) in a dataframe for further exploration. We get results for a range of k and lambdas.
## ---- to do ----
#   CHECK IF COORDINATES IMPLY WHAT SHOULD BE INSTEAD OF NA - AND THEN WHAT IS THE MEANING OF REMOVING LINKS?
# bootstrapping
# add aggregation.
# ---- results so far -----
## ---- load libraries ----
library(softImpute)
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)

## ---- functions ----

build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  library(dplyr)
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
emln_id <- 25
#layers_to_train <- 1
#layer_to_predict <- 7
#prop_ones_to_remove <- 0.2
#prop_zeros_to_remove <- 0.2

## ---- run ----
### ---- load matrices ----
# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
# Total number of layers
num_layers <- length(graph_list)

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    # Parameters for the current combination
    emln_id <- 25
    prop_ones_to_remove <- 0.2
    prop_zeros_to_remove <- 0.2
    
    # Load matrices
    d <- load_emln(emln_id)
    graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
    A_l <- d$extended
    # Total number of layers
    num_layers <- length(graph_list)
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    A[A > 0] <- 1  # Make binary
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    P[P > 0] <- 1  # Make binary
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    # # # crop to make toy example matrices
    #A <- A[1:4, 1:3]
    #P <- P[1:4, 1:3]
    
    ### ---- remove some links in P ----
    
    P_original <- P # save it for later
    num_1_to_remove <- floor(sum(P, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
    ones_in_P <- which(P == 1, arr.ind = TRUE)
    remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
    P[remove_indices] <- NA  # Set removed links to NA
    
    # Indices of zeros in P
    
    num_0_to_remove <- floor((nrow(P)*ncol(P)-sum(P, na.rm = T))*prop_zeros_to_remove)  # Number of links to remove
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    # Randomly select zeros to remove
    # set.seed(789)
    zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
    
    # Set the selected zeros to NA
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
    C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), NA, C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    # C[C>0] <- 1 # Make binary
    sum(is.na(C))
    # might need to convert C into a binary matrix
    
    ### ---- transfer learning with SVD ----
    # Define the grid of k and lambda values to search over
    k_values <- c(2, 5, 10)            # Adjust as needed
    lambda_values <- c(0.01, 0.05, 0.1)  # Adjust as needed
    
    # Initialize variables to store the best results
    results <- data.frame(k = integer(),
                          lambda = numeric(),
                          original_links = numeric(),
                          predicted_values = numeric())
    
    # Loop over all combinations of k and lambda
    for (k in k_values) {
      for (lambda in lambda_values) {
        # Apply softImpute
        fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
        
        # Reconstruct the matrix
        C_reconstructed <- softImpute::complete(C, fit)
        
        # Extract the reconstructed P matrix from C_reconstructed
        P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
        
        # Combine indices of removed ones and zeros
        test_indices <- rbind(
          data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = rep(1, nrow(remove_indices))),
          data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
        )
        
        # Get the row and column names of the test links
        test_rows <- rownames(P)[test_indices$row]  # These are the "node_to"
        test_cols <- colnames(P)[test_indices$col]  # These are the "node_from"
        
        # Actual labels and predictions
        original_links <- test_indices$label
        predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
        
        # Store the results with node information
        results <- rbind(results, data.frame(
          k = k,
          lambda = lambda,
          original_links = original_links,
          predicted_values = predicted_values,
          node_to = test_rows,
          node_from = test_cols
        ))
        
      }
    }
    
    
    ### ---- save results for current combination ----
    # After finishing the k/lambda loops, append the 'results' to 'combined_results'
    combined_results <- rbind(combined_results, data.frame(
      emln_id = emln_id,
      train_layer = layers_to_train,
      test_layer = layer_to_predict,
      prop_ones_removed = prop_ones_to_remove,
      amount_of_removed_1 = num_1_to_remove,
      prop_zeros_removed = prop_zeros_to_remove,
      amount_of_removed_0 = num_0_to_remove,
      results  # This includes k, lambda, original_links, predicted_values, node_to, node_from
    ))
  }
}

# View combined results
print(combined_results)

# Save the combined results dataframe to a CSV file
write.csv(combined_results, file = "combined_results_0.2_rem_values_names.csv", row.names = FALSE)
