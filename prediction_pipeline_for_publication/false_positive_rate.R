# add this to non-thresholded evaluation
roc_obj <- roc(response = original_binary,
               predictor = predicted_prob,
               quiet = TRUE, na.rm = TRUE,
               levels = c(0,1), direction = "<")

roc_obj <- roc(df_removed$original_binary, df_removed$predicted_prob_sigm,
               quiet = TRUE, na.rm = TRUE,
               levels = c(0,1), direction = "<")
plot(roc_obj)
auc_value <- auc(roc_obj)

threshold <- best_discrete_threshold
coords_df <- coords(
  roc_obj,
  x = threshold,
  input = "threshold",
  ret = c("specificity", "sensitivity")
)

FPR_value <- 1 - coords_df["specificity"]
TPR_value <- coords_df["sensitivity"]

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

roc_curve <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth = 1.2, color = "lightsteelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "salmon") +
  labs(
    subtitle = paste0("AUC = ", round(auc_value, 3)),
    x = "False positive rate (1 − specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  coord_equal() +
  theme_classic(base_size = 14) + tme

pdf(
  file   = "results/paper_figs/roc_curve.pdf",
  width  = 7,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(roc_curve)
dev.off()     # close the file



# add this to the evaluation
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
    f05_score = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    FPR = FP / (FP + TN), # added FPR
    nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
      sum((original_links   - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
    nnse = 1 / (2 - nse)
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f05_score = mean(f05_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    FPR = mean(FPR, na.rm = TRUE),
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()

head(result_summary)
summary(result_summary) # result_summary includes evaluation results across all iterations for each combination of islands 

# heatmap
# df_summary_fpr <- result_summary %>%
#   left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
#   rename(train_layer_name = name) %>%
#   left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
#   rename(test_layer_name = name)

# change the data frame to df_eval_summary. need to add fpr to the evaluation process
island_heatmap_fpr <- 
  ggplot(df_summary_fpr, aes(x = train_layer_name, y = test_layer_name, fill = FPR)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = df_summary_fpr[df_summary_fpr$train_layer == df_summary_fpr$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "white", mid = "lightsteelblue1", high = "steelblue", 
                       midpoint = 0.35, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added location", y = "Predicted location", fill = "False positive \nrate") +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_fpr)

pdf(
  file   = "results/paper_figs/island_heatmap_fpr.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(island_heatmap_fpr)
dev.off()     # close the file

