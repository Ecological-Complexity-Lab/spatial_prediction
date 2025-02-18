# ---- comparing scales ----
## ---- load libraries ----
library(dplyr)
library(tidyverse)
library(ggplot2)
## ---- load data ----
site_scale <- read_csv('result_summary_canary_with_distance.csv') # site scale
d <- read.csv('aggregated_equal_0_1_removal_60_1.csv') # prepare evaluators for island scale

d <- d %>%
  filter(removed == 1) %>% 
  filter(k == 2) %>% 
  filter(lambda == 0.1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0))

island_scale <- d %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_links == 1 & predicted_bin_sigm == 1),
    FN = sum(original_links == 1 & predicted_bin_sigm == 0),
    TN = sum(original_links == 0 & predicted_bin_sigm == 0),
    FP = sum(original_links == 0 & predicted_bin_sigm == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  write_csv('evaluation_df_all_itr_60_binary_island.csv')

# add distances


## ---- summarize predictive performance ----
# summarize max, min, average, se of f1, ba, recall for both scales

df1_labeled <- site_scale %>%
  mutate(dataset = "Set 1")

df2_labeled <- island_scale %>%
  mutate(dataset = "Set 2")

df_combined <- bind_rows(df1_labeled, df2_labeled)

df_long <- df_combined %>%
  pivot_longer(
    cols = c("f1_score", "recall", "mcc", "balanced_accuracy"),
    names_to = "metric",
    values_to = "value"
  )


ggplot(df_long, aes(x = metric, y = value, fill = dataset)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    title = "Comparison of Performance Metrics",
    x = "Metric",
    y = "Value"
  )


#plot
# Define custom colors for each group
custom_colors <- c("Diagonal" = "#1b9e77", 
                   "Off-diagonal (train > test)" = "#d95f02", 
                   "Off-diagonal (train < test)" = "#7570b3")

# Function to create notched boxplots with customizations
plot_f1_boxplot <- function(metric, y_axis_label = "F1 score") {
  ggplot(result_table, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = TRUE, alpha = 0.4, color = "black") +  # Notched, semi-transparent, black outline
    theme_minimal() +
    labs(x = "Layer comparison",
         y = y_axis_label) +  # Custom y-axis title
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black frame around the plot
    ) +
    tme +
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(
      c("Diagonal", "Off-diagonal (train > test)"),
      c("Diagonal", "Off-diagonal (train < test)"),
      c("Off-diagonal (train > test)", "Off-diagonal (train < test)")
    ))  # Pairwise t-tests with significance labels
}

# Generate boxplots for each F1 score metric
plot_f1_boxplot("f1_score")
