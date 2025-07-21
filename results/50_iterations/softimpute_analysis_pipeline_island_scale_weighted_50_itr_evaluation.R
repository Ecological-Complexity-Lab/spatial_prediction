# ---- softImpute analysis pipeline for island scale ----
# this pipeline allows us to take the predictions of the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# for a first time run, run first site scale and then island scale to get the working dfs with evaluators. they are both needed for scale somparison.
 
### this is the most up to date code ###
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)
library(PRROC)
library(emln)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)  # for get_legend()
library(corrplot)
library(patchwork)
library(vegan)
library(ggnewscale)
library(randomForest)
library(stringr)

#source("~/Documents/GitHub/softimpute/results/useful_for_plotting.R")
## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- parameters ----
emln_id <- 60

## ---- functions ----
# transformations

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

robust_sigmoid <- function(x, center = median(x), scale = mad(x)) {
  1 / (1 + exp(-(x - center) / scale))
}

normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

tanh_transform <- function(x){
  (tanh(x)+1)/2
}

clip_transform <- function(x){
  case_when(x<0~0,
            x>1~1,
            TRUE~x)
}

# Function to plot ROC curve with ggplot2
plot_roc_curve <- function(true_labels, predicted_scores) {
  # Create the ROC object and compute AUC
  roc_obj <- roc(response = true_labels, predictor = predicted_scores)
  auc_val <- auc(roc_obj)
  
  # Build a data frame from the ROC object for ggplot2
  df_roc <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  
  # Generate the ROC plot
  p <- ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # ROC curve line
    geom_abline(intercept = 0, slope = 1,                       # Diagonal line (random classifier)
                linetype = "dashed", color = "salmon") +
    labs(title = paste("ROC curve (AUC =", round(auc_val, 2), ")"),
         x = "False positive rate", y = "True positive rate") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

# Function to plot PR curve with ggplot2
plot_pr_curve <- function(true_labels, predicted_scores) {
  # Separate scores by class: positive (true label==1) and negative (true label==0)
  scores_pos <- predicted_scores[true_labels == 1]
  scores_neg <- predicted_scores[true_labels == 0]
  
  # Create the PR curve object; curve=TRUE returns the full curve data
  pr_obj <- pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = TRUE)
  
  # Calculate the positive class ratio for the random baseline line
  pos_ratio <- sum(true_labels == 1) / length(true_labels)
  
  # Convert the curve matrix to a data frame and set column names
  df_pr <- as.data.frame(pr_obj$curve)
  colnames(df_pr) <- c("recall", "precision", "threshold")
  
  # Generate the PR plot
  p <- ggplot(df_pr, aes(x = recall, y = precision)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # PR curve line
    geom_hline(yintercept = pos_ratio,                          # Baseline: random classifier performance
               linetype = "dashed", color = "salmon") +
    labs(title = paste("PR curve (AUC =", round(pr_obj$auc.integral, 2), ")"),
         x = "Recall", y = "Precision") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

# functions for diagonal and off-diagonal comparison
plot_boxplot <- function(data, metric, y_axis_label = "Balanced accuracy", 
                        stat_label_y = NULL, stat_size = 3) {
  
  # Calculate max value for y-axis and position for stat label
  max_y <- max(data[[metric]], na.rm = TRUE)
  label_y_position <- max_y * 0.98  # 95% of the maximum (a bit below the top)
  
  ggplot(data, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = FALSE, alpha = 0.4, color = "black") +
    theme_minimal() +
    labs(y = y_axis_label) +   # y-axis title is set via the function argument
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      axis.title.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    tme + 
    scale_fill_manual(values = custom_colors) +
    stat_compare_means(
      method = "t.test", label = "p.signif", hide.ns = FALSE, 
      comparisons = list(c("Diagonal", "Off-diagonals")),
      label.y = label_y_position * 1.1,  # Auto-set position
      size = stat_size
    ) +
    scale_y_continuous(limits = c(0.4, max_y * 1.1), labels = scales::number_format(accuracy = 0.1))  # Extend slightly above max
}

plot_hist <- function(data, metric, 
                      x_axis_label = "Balanced accuracy", 
                      y_axis_label = "Count") {
  ggplot(data, aes(x = .data[[metric]], fill = layer_comparison)) +
    geom_histogram(aes(y = ..count..), alpha = 0.4, color = "black", bins = 8, position = "dodge") +
    #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_minimal() +
    labs(x = x_axis_label,
         y = y_axis_label,
         fill = "Layer comparison") +
    theme(
      axis.text.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(values = custom_colors) + tme
}

# if you prefer density over counts:
plot_hist_density <- function(data, metric, 
                              x_axis_label = "Balanced accuracy", 
                              y_axis_label = "Density") {  
  ggplot(data, aes(x = .data[[metric]], fill = layer_comparison)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.4, color = "black", bins = 8, position = "dodge") +
    #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
    theme_minimal() +
    labs(x = x_axis_label,
         y = y_axis_label,
         fill = "Layer comparison") +
    theme(
      axis.text.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(values = custom_colors) + tme
}

# building matrices for calculating network size and density
# site scale
build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  
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

# combine plots
combine_plots <- function(p1, p2,
                          bottom_label = "Overall degree",
                          left_label = "Number of predicted, \nnon-observed links",
                          plot_margin = c(0.5, 0.5, 1, 0.3),
                          label_fontsize = 16,
                          label_fontface = "bold",
                          widths_subplots = c(1, 1),
                          final_widths = c(2, 0.3)) {
  
  # Load required packages
  require(ggplot2)
  require(gridExtra)
  require(grid)
  
  # Adjust individual plots
  p1_mod <- p1 +
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.margin = unit(plot_margin, "cm"))
  
  p2_mod <- p2 +
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.margin = unit(plot_margin, "cm"))
  
  # Arrange the two plots side-by-side
  combined_plots <- arrangeGrob(p1_mod, p2_mod, 
                                ncol = 2, 
                                widths = widths_subplots)
  
  # Add axis labels using arrangeGrob (the bottom and left text grobs)
  combined_with_axes <- arrangeGrob(
    combined_plots,
    bottom = textGrob(bottom_label, 
                      gp = gpar(fontsize = label_fontsize, fontface = label_fontface), 
                      vjust = -1.5),
    left   = textGrob(left_label, 
                      rot = 90, 
                      gp = gpar(fontsize = label_fontsize, fontface = label_fontface))
  )
  
  # Finally, arrange the whole thing with additional spacing if needed
  final_plot <- grid.arrange(
    combined_with_axes,
    ncol = 2,
    widths = final_widths
  )
  
  return(final_plot)
}

## ---- load data ----
#setwd("~/softimpute/results/results_net_60_weighted_50_itr")
df <- read_csv('weighted__scaled_island_net_60_50_itr.csv') # weighted, scaled

summary(df)

# convert negatives to zeros
df <- df %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

## ---- selecting optimal threshold ----

# 0) set thresholds
thresholds <- seq(0, 1, by = 0.1)

# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

# 2) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy= (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) /
      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)#,
    #.groups = "drop"
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup() 

# 3) average across emln_id/layer combos and pivot long
df_avg <- df_thresh %>%
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall,
      f1_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold,
               names_to  = "metric",
               values_to = "value")

# 4) plot
ggplot(df_avg, aes(threshold, value, color = metric)) +
  geom_line(size = 1) +
  labs(
    x     = "Probability threshold",
    y     = "Average metric",
    color = "Metric"
  ) +
  tme

# 1) pivot to wide so F1 and balanced_accuracy are columns
# we aim to find the optimal balance between ba and f1
df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# 2) discrete approx: minimize abs difference
best_discrete <- df_wide %>%
  mutate(absdiff = abs(f1_score - balanced_accuracy)) %>%
  slice_min(absdiff, n = 1)

# results:
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold

## ---- pr and roc curves ----
df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>% 
  mutate(predicted_prob_robust_sigm = robust_sigmoid(predicted_values)) %>% 
  mutate(predicted_minmax = normalize_min_max(predicted_values)) %>% 
  mutate(predicted_tanh = tanh_transform(predicted_values)) %>% 
  mutate(predicted_clip = clip_transform(predicted_values))

# weighted version 
# For the ROC curve:
plot_roc_curve(df_removed$original_links_binary, df_removed$predicted_values)

# For the PR curve:
plot_pr_curve(df_removed$original_links_binary, df_removed$predicted_values)

## ---- evaluation ----
### ---- plot distribution of predictions by true class ---- 
ggplot(df_removed, aes(x = predicted_prob_sigm, fill = factor(original_links_binary))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of predicted probabilities \nby true class",
       x = "Predicted probability", y = "Density",
       fill = "Original Link") +
  theme_minimal() + tme

df_removed %>% 
  filter(test_layer==4) %>% 
ggplot(aes(x = original_links, y = predicted_values)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "lightsteelblue", linetype = "dashed") +
  stat_cor(method = "spearman", label.x = 20, label.y = 50) +  # change method to "spearman" if needed
  labs(
    x = "Original links",
    y = "Predicted values",
    title = "Correlation between predictions and original links"
  ) +
  geom_abline (slope=1, linetype = "dashed", color="Red")+
  coord_equal()+
  theme_minimal(base_size = 14)

### ---- weighted version ----

df_removed <- df %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

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
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    mse = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse),
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
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()

#result_summary %>% write_csv('working_df_all_itr_60_weighted_scaled_island_50_itr.csv')
head(result_summary)
summary(result_summary)

### ---- distribution of evaluators ----
island_specificity <- ggplot(result_summary, aes(x = specificity)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  labs(x = "Specificity",
       y = "Count") +
  tme

island_f1 <- ggplot(result_summary, aes(x = f1_score)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::number_format(accuracy = 1)) +
  labs(x = "F1 score",
       y = "Count") +
  tme

island_ba <- ggplot(result_summary, aes(x = balanced_accuracy)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::number_format(accuracy = 1)) +
  labs(x = "Balanced accuracy",
       y = "Count") +
  tme

island_precision <- ggplot(result_summary, aes(x = precision)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.05)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::number_format(accuracy = 1)) +
labs(x = "Precision",
       y = "Count") +
  tme

island_recall <- ggplot(result_summary, aes(x = recall)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scales::number_format(accuracy = 1)) +
  labs(x = "Recall",
       y = "Count") +
  tme

island_rmse <- ggplot(result_summary, aes(x = rmse)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "RMSE",
       y = "Count") +
  tme

island_mse <- ggplot(result_summary, aes(x = mse)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "MSE",
       y = "Count") +
  tme

island_nse <- ggplot(result_summary, aes(x = nse)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "NSE",
       y = "Count") +
  tme

island_nnse <- ggplot(result_summary, aes(x = nnse)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "nNSE",
       y = "Count") +
  tme

p1 <- island_f1 +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p2 <- island_recall +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p3 <- island_ba +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p4 <- island_precision +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p5 <- island_specificity +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p6 <- island_rmse +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p7 <- island_mse +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p8 <- island_nse +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p9 <- island_nnse +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

combined_plots <- arrangeGrob(
  p1, p2, p3, p4, p5, p6, p7, p8, p9,
  ncol = 5, 
  nrow = 2
)
combined_with_axes <- arrangeGrob(
  combined_plots,
  #bottom = textGrob("F1 score", gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
  left   = textGrob("Count of instances", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

final_plot <- grid.arrange(
  combined_with_axes,
  ncol = 4,
  widths = c(2, 0.001, 0.2, 0)
)
