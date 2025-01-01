# ---- softImpute for predicting removed links in empirical networks ----
## ---- to do ----
# 1. make C binary? <- seems not to matter.
# fix bug! if we don't remove links we get far fewer links than we are supposed to
# we get different results for different runs (same parameters)
# 2. figure out zeros issue: how to only remove ones. in that case compare all of reconstructed P to the original
# 3. check on cropped example
## optimize thresholding: by min-max, sigmoid (logistic regression). i don't think that min-max is approapriate because the values are coordinates. i used logistic regression instead
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

# ---- results so far -----
# predicting a network without removing links with softimpute is perfect (unlike when combining a model) - also when predicting a layer with another layer
# finding optimal threshold with roc improves the predictive ability
## ---- load libraries ----
library(softImpute)
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)

## ---- functions ----

build_interaction_matrix <- function(data, layers_to_train) {
  
  # Step 1: Filter rows based on specified layers in 'layer_from'
  layers_to_filter <- paste0("layer_", layers_to_train)
  filtered_data <- subset(data, layer_from %in% layers_to_filter)
  
  # Step 2: Aggregate weights for identical species pairs
  library(dplyr)
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')
  
  # Step 3: Create the matrix with specific row and column species
  # Define unique species in node_from and node_to for matrix dimensions
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows
  
  # Initialize an empty matrix with rows as node_to species and columns as node_from species
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
layer_to_predict <- 1
prop_ones_to_remove <- 0.2
prop_zeros_to_remove <- 0

## ---- run ----
### ---- load matrices ----
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = T, directed = F)$layers_igraph # This includes in each layer all the nodes in the system (so could be many singletons)
length(graph_list) # Number of layers
unlist(lapply(graph_list, igraph::is_bipartite)) # each layer is bipartite
unlist(lapply(graph_list, igraph::vcount)) # number of nodes in each layer
A_l <- d$extended
A <- build_interaction_matrix(data = A_l, layers_to_train = layers_to_train)
dim(A)

# Organize the aggregated matrix
A[A>0] <- 1 # Make binary
A_plot <- pheatmap(A, cluster_rows = FALSE, cluster_cols = FALSE,
                   main = "Pollinator-Plant Interaction Matrix A",
                   fontsize = 8,
                   fontsize_row = 7,
                   fontsize_col = 7)

# Organize the layer to predict
P <- build_interaction_matrix(A_l, layer_to_predict)
dim(P)
P[P>0] <- 1 # Make binary

P_plot <- pheatmap(P, cluster_rows = FALSE, cluster_cols = FALSE,
                   main = "Pollinator-Plant Interaction Matrix P",
                   fontsize = 8,
                   fontsize_row = 7,
                   fontsize_col = 7)

# # # crop to make toy example matrices
#A <- A[1:4, 1:3]
#P <- P[1:4, 1:3]

### ---- remove some links in P ----

P_original <- P # save it for later
num_1_to_remove <- floor(sum(P, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
ones_in_P <- which(P == 1, arr.ind = TRUE)
remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
P[remove_indices] <- NA  # Set removed links to NA

P_plot_rem <- pheatmap(P, cluster_rows = FALSE, cluster_cols = FALSE,
                   main = "Pollinator-Plant Interaction Matrix P",
                   fontsize = 8,
                   fontsize_row = 7,
                   fontsize_col = 7)

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
    
    P_prob <- sigmoid(predicted_values)
    
    threshold <- 0.5  # Adjust based on your data
    # Determine optimal threshold
    #roc_obj <- roc(original_links, predicted_values)
    #optimal_threshold <- coords(roc_obj, "best", ret="threshold")
    #optimal_threshold <- as.numeric(optimal_threshold)
    #predicted_links <- ifelse(P_prob > threshold, 1, 0)
    # Classify links based on probabilities in P_prob
    predicted_links <- ifelse(P_prob == threshold, 0, 1)
    
    
    # Confusion matrix
    confusion <- table(Predicted = predicted_links, Actual = original_links)
    print(confusion)
    
    # # Accuracy
    # accuracy <- sum(predicted_links == original_links) / length(predicted_links)
    # cat("Accuracy:", accuracy, "\n")
    
    # Precision, Recall, F1 Score
    # Confusion matrix components
    # Ensure both vectors have no NA values (since P_original should have no NAs)
    # If P_original still contains NAs, remove them or handle appropriately
    
    TP <- sum(predicted_links == 1 & original_links == 1)  # True Positives
    FP <- sum(predicted_links == 1 & original_links == 0)  # False Positives
    TN <- sum(predicted_links == 0 & original_links == 0)  # True Negatives
    FN <- sum(predicted_links == 0 & original_links == 1)  # False Negatives

    
    # Calculate Sensitivity (Recall) and Specificity
    Sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    Specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
    
    # Calculate Balanced Accuracy
    Balanced_Accuracy <- (Sensitivity + Specificity) / 2
    
    # Additional metrics
    Precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
    Recall <- Sensitivity
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

# Output the best parameters and corresponding balanced accuracy
cat("Best Balanced Accuracy:", best_balanced_accuracy, "\n")
cat("Best k (rank.max):", best_k, "\n")
cat("Best lambda:", best_lambda, "\n")

results$num_1_to_remove <- num_1_to_remove
results$num_0_to_remove <- num_0_to_remove

# View all results
print(results)

# 
df <- data.frame(
  emln_id = emln_id,  # Example IDs
  layers_to_train = layers_to_train,  # Layers to train
  layer_to_predict = layer_to_predict,  # Layers to predict
  prop_ones_to_remove = prop_ones_to_remove,
  TP = TP,
  TN = TN,
  FP = FP,
  FN = FN,
  specificity = Specificity,
  precision = Precision,
  recall = Recall,
  f1_score = F1_Score,
  balanced_accuracy = best_balanced_accuracy,
  prop_zeros_to_remove = prop_zeros_to_remove
)

view(df)
# Optional: Plot Balanced Accuracy as a function of k and lambda
ggplot(results, aes(x = factor(k), y = Balanced_Accuracy, color = factor(lambda))) +
  geom_point(size = 3) +
  geom_line(aes(group = factor(lambda))) +
  labs(title = "Balanced Accuracy for different k and lambda values",
       x = "k (rank.max)",
       y = "Balanced Accuracy",
       color = "lambda") +
  theme_minimal()

# Show the row with the highest Balanced Accuracy
results[results$Balanced_Accuracy == max(results$Balanced_Accuracy),]


### ---- try to predict P ----
# lambda=0.01
# k=8
fit <- softImpute(P, rank.max = best_k, lambda = best_lambda, type = "svd", )
P_reconstructed <- softImpute::complete(P, fit)
# P_reconstructed[P_reconstructed>0] <- 1 # Make binary
# P_reconstructed[P_reconstructed<0] <- 0 # Make binary # convert to binary after solving the threshold issue

P_reconstructed <- sigmoid(P_reconstructed)
P_reconstructed <- ifelse(P_reconstructed == threshold, 0, 1)

P_predicted <- pheatmap(P_reconstructed, cluster_rows = FALSE, cluster_cols = FALSE,
                        main = "Pollinator-Plant Interaction Matrix P (predicted)",
                        fontsize = 8,
                        fontsize_row = 7,
                        fontsize_col = 7)

# Combine the plots side by side
grid.arrange(P_predicted$gtable, P_plot_rem$gtable, A_plot$gtable, 
             ncol = 3)
