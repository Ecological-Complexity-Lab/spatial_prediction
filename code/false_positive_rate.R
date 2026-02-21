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

ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth = 1.2, color = "lightsteelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "salmon") +
  labs(
    title = "ROC curve",
    subtitle = paste0("AUC = ", round(auc_value, 3)),
    x = "False positive rate (1 − specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  coord_equal() +
  theme_classic(base_size = 14) + tme
