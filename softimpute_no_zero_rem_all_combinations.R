# ---- softImpute for predicting removed links in empirical networks - all layer combinations ----
## ---- to do ----
#   CHECK IF COORDINATES IMPLY WHAT SHOULD BE INSTEAD OF NA - AND THEN WHAT IS THE MEANING OF REMOVING LINKS?
# compare the links that i did not remove in the reconstructed matrix to the original ones. if they are perfectly predicted it's not prediction
# 1. make C binary? <- seems not to matter.
# fix bug! if we don't remove links we get far fewer links than we are supposed to
# we get different results for different runs (same parameters)
# 2. figure out zeros issue: how to only remove ones. in that case compare all of reconstructed P to the original
# 3. check on cropped example
## optimize thresholding: by min-max, sigmoid (logistic regression). i don't think that min-max is appropriate because the values are coordinates. i used logistic regression instead
# C reconstructed when A=P should be identical to P and A, check <- it's okay
# 4. figure out whats lambda
# 5. try to maximize recall instead of balanced accuracy
# 6. fix P_reconstructed values. some of them are negative <- i started by using roc analysis to find a more appropriate threshold.
# 7. check evaluation according to the goal
# 8. why is BA decrease with increasing k? in the toy example it increases
# 9. check proportions of shared species
# 10. wrap in a function and analyze for each pair of layers. <- done in softimpute_empiric_all_combinations. should be adjusted according to the changes in this code.
# 11. bootstrapping
# according to chat: The softImpute method is designed for continuous data and may not be optimal for binary matrices. so i tried to apply on non-binary matrices, which resulted in decreased predictive ability.
# problem. removed links that are not in A are still predicted.
# ---- results so far -----
# F1 score is significantly better for the diagonal than when adding information from other layers, while the recall of removed links is significantly better when incorporating another layer. which is strange because the recalled links are sometimes not even present in the added layer. balanced accuracy is similar. adding information from other layers might lower overfitting. When A=P, the model may overfit to the patterns of P, resulting in high precision but possibly missing out on discovering new links. When A is different than P, the model generalizes better by incorporating diverse patterns but risks introducing errors due to discrepancies between the networks.
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
layers_to_train <- 1
layer_to_predict <- 7
prop_ones_to_remove <- 0.2
prop_zeros_to_remove <- 0

## ---- run ----
### ---- load matrices ----
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
    prop_zeros_to_remove <- 0
    
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
    
    pheatmap(C, cluster_rows = FALSE, cluster_cols = FALSE,
             main = "Pollinator-Plant Interaction Matrix C")
    
    ### ---- transfer learning with SVD ----
    # Define the grid of k and lambda values to search over
    k_values <- seq(2, 10, by = 2)            # Adjust as needed
    lambda_values <- seq(0.01, 0.1, by = 0.03)  # Adjust as needed
    
    # Initialize variables to store the best results
    best_balanced_accuracy <- 0
    best_k <- NULL
    best_lambda <- NULL
    results <- data.frame(k = integer(),
                          lambda = numeric(),
                          Balanced_Accuracy = numeric(),
                          Sensitivity = numeric(),
                          Specificity = numeric(),
                          Precision = numeric(),
                          Recall = numeric(),
                          F1_Score = numeric())
    
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
        
        # Get the row and column names
        test_rows <- rownames(P)[test_indices$row]
        test_cols <- colnames(P)[test_indices$col]
        
        # Actual labels
        original_links <- test_indices$label
        predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
        
        # Threshold to decide if a link is present
        # first convert predicted values to probabilities
        # Apply sigmoid transformation
        sigmoid <- function(x) {
          1 / (1 + exp(-x))
        }
        
        P_prob <- sigmoid(P_reconstructed)
        
        # for removed links only
        P_prob_rem <- sigmoid(predicted_values)
        
        threshold <- 0.5  # Adjust based on your data
        # Determine optimal threshold
        #roc_obj <- roc(original_links, predicted_values)
        #optimal_threshold <- coords(roc_obj, "best", ret="threshold")
        #optimal_threshold <- as.numeric(optimal_threshold)
        #predicted_links <- ifelse(P_prob > threshold, 1, 0)
        # Classify links based on probabilities in P_prob
        predicted_links <- ifelse(P_prob == threshold, 0, 1)
        predicted_rem_links <- ifelse(P_prob_rem == threshold, 0, 1)
        
        # Confusion matrix
        confusion <- table(Predicted = predicted_links, Actual = P_original)
        print(confusion)
        
        # # Accuracy
        # accuracy <- sum(predicted_links == original_links) / length(predicted_links)
        # cat("Accuracy:", accuracy, "\n")
        
        # Precision, Recall, F1 Score
        # Confusion matrix components
        # Ensure both vectors have no NA values (since P_original should have no NAs)
        # If P_original still contains NAs, remove them or handle appropriately
        
        # Compute confusion matrix components
        TP <- sum(P_original == 1 & predicted_links == 1)
        TN <- sum(P_original == 0 & predicted_links == 0)
        FP <- sum(P_original == 0 & predicted_links == 1)
        FN <- sum(P_original == 1 & predicted_links == 0)
        
        # for removed links specifically:
        TP_rem <- sum(predicted_rem_links == 1 & original_links == 1)  # True Positives
        FP_rem <- sum(predicted_rem_links == 1 & original_links == 0)  # False Positives
        TN_rem <- sum(predicted_rem_links == 0 & original_links == 0)  # True Negatives
        FN_rem <- sum(predicted_rem_links == 0 & original_links == 1)  # False Negatives
        
        # Calculate Sensitivity (Recall) and Specificity
        Sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
        Specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
        
        # Calculate Balanced Accuracy
        Balanced_Accuracy <- (Sensitivity + Specificity) / 2
        
        # Additional metrics
        Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
        Recall <- Sensitivity
        Recall_rem <- ifelse((TP_rem + FN_rem) > 0, TP_rem / (TP_rem + FN_rem), 0)
        F1_Score <- ifelse((Precision + Recall) > 0, 2 * (Precision * Recall) / (Precision + Recall), 0)
        
        # Store the results
        results <- rbind(results, data.frame(k = k,
                                             lambda = lambda,
                                             Balanced_Accuracy = Balanced_Accuracy,
                                             Sensitivity = Sensitivity,
                                             Specificity = Specificity,
                                             Precision = Precision,
                                             Recall = Recall,
                                             F1_Score = F1_Score))
        
        # Update the best parameters if current balanced accuracy is higher
        if (Balanced_Accuracy > best_balanced_accuracy) {
          best_balanced_accuracy <- Balanced_Accuracy
          best_k <- k
          best_lambda <- lambda
        }
      }
    }
    
    
    ### ---- evaluation ----
    # Output the best parameters and corresponding balanced accuracy
    cat("Best Balanced Accuracy:", best_balanced_accuracy, "\n")
    cat("Best k (rank.max):", best_k, "\n")
    cat("Best lambda:", best_lambda, "\n")
    
    ### ---- save results for current combination ----
    combined_results <- rbind(combined_results, data.frame(
      emln_id = emln_id,  # Example IDs
      layers_to_train = layers_to_train,  # Layers to train
      layer_to_predict = layer_to_predict,  # Layers to predict
      prop_ones_to_remove = prop_ones_to_remove,
      amount_of_removed_1 = num_1_to_remove,
      prop_zeros_to_remove = prop_zeros_to_remove,
      amount_of_removed_0 = num_0_to_remove,
      TP = TP,
      TN = TN,
      FP = FP,
      FN = FN,
      specificity = Specificity,
      precision = Precision,
      recall = Recall,
      recall_rem = Recall_rem,
      f1_score = F1_Score,
      balanced_accuracy = best_balanced_accuracy
    ))
    }
}

# View combined results
print(combined_results)

# Save the combined results dataframe to a CSV file
write.csv(combined_results, file = "combined_results_no_0_rem.csv", row.names = FALSE)

## ---- visualize ----
# produce a heatmap
combined_results$diagonal <- combined_results$layers_to_train == combined_results$layer_to_predict

layer_to_layer_plot <- ggplot(combined_results, aes(x = layers_to_train, y = layer_to_predict, fill = balanced_accuracy)) +
  # Entire heatmap with white borders
  geom_tile(color = "white", size = 0.1) +  
  # Highlight diagonal tiles with black borders
  geom_tile(data = combined_results[combined_results$diagonal, ],
            aes(x = layers_to_train, y = layer_to_predict, fill = balanced_accuracy), 
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
    panel.grid.minor = element_blank()) +      # Remove minor grid lines
    #axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  # Rotate x-axis labels
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure each x tick has a label
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  coord_fixed()  # Ensure tiles are square

# Display the plot
print(layer_to_layer_plot)

## ---- test differences ----
results <- read.csv("combined_results_no_0_rem.csv")

# Separate data into two groups
same_layers <- results$balanced_accuracy[results$layers_to_train == results$layer_to_predict]
# Selecting cases where layer_to_train is different from layer_to_predict
different_layers <- results$balanced_accuracy[results$layers_to_train != results$layer_to_predict]

# results$layers_match <- if(results$layers_to_train == results$layer_to_predict
#                            results$layers_match == TRUE
#                            ifelse(results$layers_match == FALSE)

# Add a new column indicating whether the layers match
results$layers_match <- results$layers_to_train == results$layer_to_predict

# Boxplot for visual comparison

boxplot(balanced_accuracy ~ layers_match, data = results,
        #main = "Comparison of Balanced Accuracy",
        xlab = "Layers match (same vs different)",
        ylab = "Balanced accuracy",
        col = c("skyblue", "lightgreen"))

# Perform a two-sample t-test
t_test <- t.test(same_layers, different_layers, var.equal = FALSE)
print(t_test)
