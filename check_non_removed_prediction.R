
    
    # Identify non-removed indices
    non_removed_indices <- which(!is.na(P), arr.ind = TRUE)
    
    # Actual values from P_original at non-removed positions
    actual_non_removed <- P_original[non_removed_indices]
    
    # Predicted values from P_reconstructed at the same positions
    predicted_non_removed <- P_reconstructed[non_removed_indices]
    
    identical(actual_non_removed, predicted_non_removed)
    
    # Identify discrepancies
    discrepancies <- which(actual_non_removed != predicted_non_removed)
    num_discrepancies <- length(discrepancies)
    total_non_removed <- length(actual_non_removed)
    discrepancy_rate <- num_discrepancies / total_non_removed
    
    # Get the indices of discrepancies in terms of rows and columns
    discrepancy_indices <- non_removed_indices[discrepancies, , drop = FALSE]
    
    # Get the species names (row and column names)
    discrepancy_rows <- rownames(P)[discrepancy_indices[, "row"]]
    discrepancy_cols <- colnames(P)[discrepancy_indices[, "col"]]
    
    # Create a data frame with details of discrepancies
    discrepancy_details <- data.frame(
      row_index = discrepancy_indices[, "row"],
      col_index = discrepancy_indices[, "col"],
      species_row = discrepancy_rows,
      species_col = discrepancy_cols,
      actual = actual_non_removed[discrepancies],
      predicted = predicted_non_removed[discrepancies]
    )
    
