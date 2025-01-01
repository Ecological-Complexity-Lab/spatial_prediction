# Install and load the softImpute package
# install.packages("softImpute")
library(softImpute)
library(PRROC)
library(pROC)
library(ggplot2)

# Set seed for reproducibility
# set.seed(123)

# Dimensions of matrices A and P
m_A <- 20  # Number of rows in A
n_A <- 40  # Number of columns in A
m_P <- 15  # Number of rows in P
n_P <- 35  # Number of columns in P

# Number of shared nodes
shared_rows <- 8   # Shared rows between A and P
shared_cols <- 12  # Shared columns between A and P

# Generate unique identifiers for nodes
A_rows <- paste0("A_row_", 1:(m_A - shared_rows))
A_cols <- paste0("A_col_", 1:(n_A - shared_cols))
P_rows <- paste0("P_row_", 1:(m_P - shared_rows))
P_cols <- paste0("P_col_", 1:(n_P - shared_cols))
shared_rows_ids <- paste0("Shared_row_", 1:shared_rows)
shared_cols_ids <- paste0("Shared_col_", 1:shared_cols)

# Combine node identifiers
A_row_ids <- c(A_rows, shared_rows_ids)
A_col_ids <- c(A_cols, shared_cols_ids)
P_row_ids <- c(P_rows, shared_rows_ids)
P_col_ids <- c(P_cols, shared_cols_ids)

# Create empty matrices A and P
A <- matrix(0, nrow = length(A_row_ids), ncol = length(A_col_ids),
            dimnames = list(A_row_ids, A_col_ids))
P <- matrix(0, nrow = length(P_row_ids), ncol = length(P_col_ids),
            dimnames = list(P_row_ids, P_col_ids))

# Function to generate binary data with higher probability for shared nodes
generate_data <- function(rows, cols, shared_rows, shared_cols, prob_shared, prob_non_shared) {
  data_matrix <- matrix(0, nrow = length(rows), ncol = length(cols),
                        dimnames = list(rows, cols))
  
  for (i in rows) {
    for (j in cols) {
      if (i %in% shared_rows && j %in% shared_cols) {
        # Shared nodes: higher probability to be the same
        data_matrix[i, j] <- rbinom(1, 1, prob_shared)
      } else {
        # Non-shared nodes
        data_matrix[i, j] <- rbinom(1, 1, prob_non_shared)
      }
    }
  }
  return(data_matrix)
}

# Generate data for A and P
A <- generate_data(A_row_ids, A_col_ids, shared_rows_ids, shared_cols_ids, prob_shared = 0.6, prob_non_shared = 0.1)
P <- generate_data(P_row_ids, P_col_ids, shared_rows_ids, shared_cols_ids, prob_shared = 0.6, prob_non_shared = 0.1)

# Remove some links in P (simulate missing links)
# set.seed(456)
prop_ones_to_remove <- 0.2
num_1_to_remove <- floor(sum(P, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
ones_in_P <- which(P == 1, arr.ind = TRUE)
remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
P[remove_indices] <- NA  # Set removed links to NA


# Indices of zeros in P
prop_zeros_to_remove <- 0.01
num_0_to_remove <- floor((nrow(P)*ncol(P)-sum(P, na.rm = T))*prop_zeros_to_remove)  # Number of links to remove
zeros_in_P <- which(P == 0, arr.ind = TRUE)
# Randomly select zeros to remove
# set.seed(789)
zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]

# Set the selected zeros to NA
P[zeros_to_remove_indices] <- NA


# Combine A and P into a single matrix C with NAs representing missing data
all_row_ids <- unique(c(A_row_ids, P_row_ids))
all_col_ids <- unique(c(A_col_ids, P_col_ids))
C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
            dimnames = list(all_row_ids, all_col_ids))

# Place A into C
C[A_row_ids, A_col_ids] <- A

# Place P into C
# Ensure that existing entries are not overwritten; sum overlapping entries
C[P_row_ids, P_col_ids] <- ifelse(is.na(C[P_row_ids, P_col_ids]), NA, C[P_row_ids, P_col_ids] + P[P_row_ids, P_col_ids])
sum(is.na(C))


# Transfer Learning with SVD on Toy Examples
# Your data generation code remains the same
# [Your existing code for data generation and setting up matrices A, P, and C]

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
    fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 500)
    
    # Reconstruct the matrix
    C_reconstructed <- softImpute::complete(C, fit)
    
    # Extract the reconstructed P matrix from C_reconstructed
    P_reconstructed <- C_reconstructed[P_row_ids, P_col_ids]
    
    
    # Combine indices of removed ones and zeros
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = 0)
    )
    
    # Get the row and column names
    test_rows <- rownames(P)[test_indices$row]
    test_cols <- colnames(P)[test_indices$col]
    
    # Actual labels
    original_links <- test_indices$label
    predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
    
    # Threshold to decide if a link is present
    threshold <- 0.5  # Adjust based on your data
    predicted_links <- ifelse(predicted_values > threshold, 1, 0)
    # Confusion matrix
    confusion <- table(Predicted = predicted_links, Actual = original_links)
    print(confusion)
    
    # # Accuracy
    # accuracy <- sum(predicted_links == original_links) / length(predicted_links)
    # cat("Accuracy:", accuracy, "\n")
    
    # Precision, Recall, F1 Score
    # Confusion matrix components
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


# Try to predict P
lambda=0.01
k=8
fit <- softImpute(P, rank.max = k, lambda = lambda, type = "svd", )
P_reconstructed <- softImpute::complete(P, fit)

