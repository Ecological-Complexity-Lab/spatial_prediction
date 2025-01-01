# ---- softImpute for predicting removed links in empirical networks: for each combination of layers ----
## ---- to do ----
# 1. make C binary? <- seems not to matter.
# 2. figure out zeros issue: how to only remove ones. in that case compare all of reconstructed P to the original
# 3. check on cropped example
# figure out whats lambda
# try to maximize recall instead of balanced accuracy
# fix P_reconstructed values. some of them are negative
# 4. check evaluation according to the goal
# why is BA decrease with increasing k? in the toy example it increases
# check proportions of shared species
# wrap in a function and analyze for each pair of layers.
# how to infer labels from predicted values??
# 5. bootstrapping

## ---- load libraries ----
# install.packages("softImpute")
library(softImpute)
library(PRROC)
library(pROC)
library(ggplot2)
library(emln)
library(pheatmap)

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


# ## ---- parameters ----
# emln_id <- 25
# prop_ones_to_remove <- 0.2
# prop_zeros_to_remove <- 0.01

## ---- run ----

# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

# Total number of layers
num_layers <- length(graph_list)

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    # Parameters for the current combination
    emln_id <- 25
    prop_ones_to_remove <- 0.2
    prop_zeros_to_remove <- 0.01
    
    # Load matrices
    d <- load_emln(emln_id)
    graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
    A_l <- d$extended
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    A[A > 0] <- 1  # Make binary
    
    # Build the layer to predict matrix P
    P <- P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    P[P > 0] <- 1  # Make binary
    
    ### ---- Remove some links in P ----
    num_1_to_remove <- floor(sum(P, na.rm = TRUE) * prop_ones_to_remove)
    ones_in_P <- which(P == 1, arr.ind = TRUE)
    remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
    P[remove_indices] <- NA  # Set removed links to NA
    
    num_0_to_remove <- floor((nrow(P) * ncol(P) - sum(P, na.rm = TRUE)) * prop_zeros_to_remove)
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
    P[zeros_to_remove_indices] <- NA  # Set zeros to NA
    
    ### ---- Create combined matrix C ----
    all_row_ids <- unique(c(rownames(A), rownames(P)))
    all_col_ids <- unique(c(colnames(A), colnames(P)))
    C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                dimnames = list(all_row_ids, all_col_ids))
    C[rownames(A), colnames(A)] <- A
    C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), NA, C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    
    ### ---- Transfer learning with SVD ----
    k_values <- seq(2, 10, by = 2)
    lambda_values <- seq(0.01, 0.1, by = 0.03)
    best_balanced_accuracy <- 0
    results <- data.frame()
    
    for (k in k_values) {
      for (lambda in lambda_values) {
        fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
        C_reconstructed <- softImpute::complete(C, fit)
        P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
        
        # Combine indices of removed ones and zeros
        test_indices <- rbind(
          data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = 1),
          data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = 0)
        )
        
        test_rows <- rownames(P)[test_indices$row]
        test_cols <- colnames(P)[test_indices$col]
        original_links <- test_indices$label
        predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
        predicted_links <- ifelse(predicted_values > 0.5, 1, 0)
        
        TP <- sum(predicted_links == 1 & original_links == 1)
        FP <- sum(predicted_links == 1 & original_links == 0)
        TN <- sum(predicted_links == 0 & original_links == 0)
        FN <- sum(predicted_links == 0 & original_links == 1)
        Sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
        Specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
        Balanced_Accuracy <- (Sensitivity + Specificity) / 2
        Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
        Recall <- Sensitivity
        F1_Score <- ifelse((Precision + Recall) > 0, 2 * (Precision * Recall) / (Precision + Recall), 0)
        
        results <- rbind(results, data.frame(
          k = k, lambda = lambda, Balanced_Accuracy = Balanced_Accuracy,
          Sensitivity = Sensitivity, Specificity = Specificity,
          Precision = Precision, Recall = Recall, F1_Score = F1_Score
        ))
        
        if (Balanced_Accuracy > best_balanced_accuracy) {
          best_balanced_accuracy <- Balanced_Accuracy
        }
      }
    }
    
    ### ---- Save Results for Current Combination ----
    combined_results <- rbind(combined_results, data.frame(
      emln_id = emln_id,
      layers_to_train = layers_to_train,
      layer_to_predict = layer_to_predict,
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      specificity = Specificity,
      precision = Precision,
      recall = Recall,
      f1_score = F1_Score,
      best_balanced_accuracy = best_balanced_accuracy,
      prop_ones_to_remove = prop_ones_to_remove,
      prop_zeros_to_remove = prop_zeros_to_remove
    ))
  }
}

# View combined results
print(combined_results)

# Save the combined results dataframe to a CSV file
write.csv(combined_results, file = "combined_results.csv", row.names = FALSE)

### ---- test differences ----
results <- read.csv("combined_results_softimpute.csv")

# Separate data into two groups
same_layers <- results$best_balanced_accuracy[results$layers_to_train == results$layer_to_predict]
# Selecting cases where layer_to_train is different from layer_to_predict
different_layers <- results$best_balanced_accuracy[results$layers_to_train != results$layer_to_predict]

# results$layers_match <- if(results$layers_to_train == results$layer_to_predict
#                            results$layers_match == TRUE
#                            ifelse(results$layers_match == FALSE)

                           # Add a new column indicating whether the layers match
results$layers_match <- results$layers_to_train == results$layer_to_predict
                           
# Boxplot for visual comparison

boxplot(best_balanced_accuracy ~ layers_match, data = results,
        #main = "Comparison of Balanced Accuracy",
        xlab = "Layers match (same vs different)",
        ylab = "Balanced accuracy",
        col = c("skyblue", "lightgreen"))

# Perform a two-sample t-test
t_test <- t.test(same_layers, different_layers, var.equal = FALSE)
print(t_test)

# produce a heatmap
results$diagonal <- results$layers_to_train == results$layer_to_predict

layer_to_layer_plot <- ggplot(results, aes(x = layers_to_train, y = layer_to_predict, fill = best_balanced_accuracy)) +
  # Entire heatmap with white borders
  geom_tile(color = "white", size = 0.1) +  
  # Highlight diagonal tiles with black borders
  geom_tile(data = results[results$diagonal, ],
            aes(x = layers_to_train, y = layer_to_predict, fill = best_balanced_accuracy), 
            color = "black", size = 1.2) +  
  # Gradient color scale for balanced accuracy
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  
  # Add axis labels and title for color bar
  labs(x = "Training layer", y = "Predicted layer", fill = "Balanced\naccuracy") +
  # Minimal theme for better appearance
  theme_minimal() +
  # Customize theme elements
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Reduce plot margins
    panel.background = element_blank(),       # No panel background
    panel.grid.major = element_blank(),       # Remove major grid lines
    panel.grid.minor = element_blank(),       # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  # Rotate x-axis labels
  #   axis.title.x = element_text(linewidth = 12),  # Adjust x-axis title size
  #   axis.title.y = element_text(linewidth = 12),  # Adjust y-axis title size
  #   legend.title = element_text(linewidth = 10),  # Adjust legend title size
  #   legend.text = element_text(linewidth = 8)     # Adjust legend text size
  # ) +
  coord_fixed()  # Ensure tiles are square

# Display the plot
print(layer_to_layer_plot)


