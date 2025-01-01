library(pROC)

# Initialize variables to store the best ROC-based results
best_balanced_accuracy <- 0
best_roc_auc <- 0
best_roc_k <- NULL
best_roc_lambda <- NULL
best_roc_threshold <- NULL
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
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = 0)
    )
    
    # Get the row and column names
    test_rows <- rownames(P)[test_indices$row]
    test_cols <- colnames(P)[test_indices$col]
    
    # Actual labels
    original_links <- test_indices$label
    predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]  # Predicted scores as a vector
    
    # Ensure lengths match
    if (length(predicted_values) != length(original_links)) {
      stop("Mismatch in lengths of predicted and original labels!")
    }
    
    # Perform ROC analysis
    roc_curve <- roc(original_links, predicted_values)
    optimal_threshold <- coords(roc_curve, "best", ret = "threshold")
    auc_value <- auc(roc_curve)
    
    # Threshold-based classification
    predicted_links <- ifelse(predicted_values > optimal_threshold, 1, 0)
    
    # Confusion matrix
    #confusion <- table(Predicted = predicted_links, Actual = original_links)
    #print(confusion)
    
    # Calculate metrics based on the optimal threshold
    TP <- sum(predicted_links == 1 & original_links == 1)  # True Positives
    FP <- sum(predicted_links == 1 & original_links == 0)  # False Positives
    TN <- sum(predicted_links == 0 & original_links == 0)  # True Negatives
    FN <- sum(predicted_links == 0 & original_links == 1)  # False Negatives
    
    Sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
    Specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
    Balanced_Accuracy <- (Sensitivity + Specificity) / 2
    
    # Store the ROC-based results
    if (auc_value > best_roc_auc) {
      best_roc_auc <- auc_value
      best_roc_k <- k
      best_roc_lambda <- lambda
      best_roc_threshold <- optimal_threshold
    }
    
    cat("AUC:", auc_value, "Optimal Threshold:", optimal_threshold, "\n")
  }
}

# Output the best ROC-based parameters and corresponding AUC
cat("Best AUC:", best_roc_auc, "\n")
cat("Best k (rank.max):", best_roc_k, "\n")
cat("Best lambda:", best_roc_lambda, "\n")
cat("Best ROC Threshold:", best_roc_threshold, "\n")
