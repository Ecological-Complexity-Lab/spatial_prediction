layer_to_layer_plot_f1 <- ggplot(combined_results, aes(x = layers_to_train, y = layer_to_predict, fill = f1_score)) +
  # Entire heatmap with white borders
  geom_tile(color = "white", size = 0.1) +  
  # Highlight diagonal tiles with black borders
  geom_tile(data = combined_results[combined_results$diagonal, ],
            aes(x = layers_to_train, y = layer_to_predict, fill = f1_score), 
            color = "black", size = 1.2) +  
  # Gradient color scale for balanced accuracy
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  
  # Add axis labels and title for color bar
  labs(x = "Training layer", y = "Predicted layer", fill = "F1\nscore") +
  # Minimal theme for better appearance
  theme_minimal() +
  # Customize theme elements
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Reduce plot margins
    panel.background = element_blank(),       # No panel background
    panel.grid.major = element_blank(),       # Remove major grid lines
    panel.grid.minor = element_blank(),       # Remove minor grid lines
    axis.text.x = element_text(hjust = 0.5, vjust = 1)) +  # Rotate x-axis labels
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure each x tick has a label
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  coord_fixed()  # Ensure tiles are square

# Display the plot
print(layer_to_layer_plot_f1)

## ---- test differences ----
results <- read.csv("combined_results_no_0_rem.csv")

# Separate data into two groups
same_layers <- results$f1_score[results$layers_to_train == results$layer_to_predict]
# Selecting cases where layer_to_train is different from layer_to_predict
different_layers <- results$f1_score[results$layers_to_train != results$layer_to_predict]

# results$layers_match <- if(results$layers_to_train == results$layer_to_predict
#                            results$layers_match == TRUE
#                            ifelse(results$layers_match == FALSE)

# Add a new column indicating whether the layers match
results$layers_match <- results$layers_to_train == results$layer_to_predict

# Boxplot for visual comparison

boxplot(f1_score ~ layers_match, data = results,
        #main = "Comparison of Balanced Accuracy",
        xlab = "",
        ylab = "F1 score",
        col = c("skyblue", "lightgreen"))

# Perform a two-sample t-test
t_test <- t.test(same_layers, different_layers, var.equal = FALSE)
print(t_test)

layer_to_layer_plot_recall_rem <- ggplot(combined_results, aes(x = layers_to_train, y = layer_to_predict, fill = recall_rem)) +
  # Entire heatmap with white borders
  geom_tile(color = "white", size = 0.1) +  
  # Highlight diagonal tiles with black borders
  geom_tile(data = combined_results[combined_results$diagonal, ],
            aes(x = layers_to_train, y = layer_to_predict, fill = recall_rem), 
            color = "black", size = 1.2) +  
  # Gradient color scale for balanced accuracy
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  
  # Add axis labels and title for color bar
  labs(x = "Training layer", y = "Predicted layer", fill = "recall\nremoved") +
  # Minimal theme for better appearance
  theme_minimal() +
  # Customize theme elements
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Reduce plot margins
    panel.background = element_blank(),       # No panel background
    panel.grid.major = element_blank(),       # Remove major grid lines
    panel.grid.minor = element_blank(),       # Remove minor grid lines
    axis.text.x = element_text(hjust = 0.5, vjust = 1)) +  # Rotate x-axis labels
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure each x tick has a label
  scale_y_continuous(breaks = seq(1, 7, by = 1)) +
  coord_fixed()  # Ensure tiles are square

# Display the plot
print(layer_to_layer_plot_recall_rem)

## ---- test differences ----
results <- read.csv("combined_results_no_0_rem.csv")

# Separate data into two groups
same_layers <- results$recall_rem[results$layers_to_train == results$layer_to_predict]
# Selecting cases where layer_to_train is different from layer_to_predict
different_layers <- results$recall_rem[results$layers_to_train != results$layer_to_predict]

# Add a new column indicating whether the layers match
results$layers_match <- results$layers_to_train == results$layer_to_predict

# Boxplot for visual comparison

boxplot(recall_rem ~ layers_match, data = results,
        xlab = "predicted = train",
        ylab = "recall of removed links",
        col = c("skyblue", "lightgreen"))

# Perform a two-sample t-test
t_test <- t.test(same_layers, different_layers, var.equal = FALSE)
print(t_test)

