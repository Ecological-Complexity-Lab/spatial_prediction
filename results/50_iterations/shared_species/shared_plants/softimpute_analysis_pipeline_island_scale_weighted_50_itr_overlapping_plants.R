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
#setwd("~/Documents/github/softimpute/results/50_iterations/shared_species")
df <- read_csv('weighted__scaled_island_net_60_50_itr_shared_plants.csv') # weighted, scaled

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

combined_plots <- arrangeGrob(
  p1, p2, p3, p4, p5, p6, p7,
  ncol = 4, 
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

## ---- how many iterations we had? ----
# Count number of iterations for each combination
itr_counts <- df_removed %>%
  group_by(train_layer, test_layer) %>%
  summarise(n_itr = n_distinct(itr)) %>%
  ungroup()

# Plot it
ggplot(itr_counts, aes(x = factor(train_layer), y = n_itr, fill = factor(test_layer))) +
  geom_col(position = "dodge") +
  labs(
    title = "Number of Iterations per Train/Test Layer",
    x = "Train Layer",
    y = "Number of Iterations",
    fill = "Test Layer"
  ) +
  theme_minimal()


## ---- diagonal vs. off-diagonals ----
# this analysis shows us if predictions made using added information from other locations (off-diagonals in layer-to-layer predictions, as a heatmap) is any better than not adding any information (cases on the diagonal)

result_summary <- result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Diagonal",
    train_layer != test_layer ~ "Off-diagonals"
  ))

custom_colors <- c("Diagonal" = "steelblue",
                   "Off-diagonals" = "thistle")

# boxplots:
island_ba <- plot_boxplot(result_summary, metric = "balanced_accuracy", 
                        y_axis_label = "Balanced accuracy", stat_label_y = 0.85, stat_size = 6)
island_f1 <- plot_boxplot(result_summary, metric = "f1_score", 
                        y_axis_label = "F1 score", stat_label_y = 0.7, stat_size = 5)

island_mse_b <- plot_boxplot(result_summary, metric = "mse", 
                          y_axis_label = "MSE", stat_size = 6)

combined_plots <- arrangeGrob(
  island_ba, island_f1, 
  ncol = 2, 
  widths = c(1, 1)
)

final_plot <- grid.arrange(
  combined_plots,
  ncol = 2,
  widths = c(2, 0.3)
)

plot_boxplot <- function(data, metric, y_axis_label = "Balanced accuracy", 
                         stat_size = 3) {
  
  # Calculate max value for y-axis and position for stat label
  max_y <- max(data[[metric]], na.rm = TRUE)
  min_y <- min(data[[metric]], na.rm = TRUE)
  label_y_position <- max_y * 0.92  # 95% of the maximum (a bit below the top)
  
  ggplot(data, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = FALSE, alpha = 0.4, color = "black") +
    theme_minimal() +
    labs(y = y_axis_label) +
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
    scale_y_continuous(limits = c(min_y, max_y * 1.2), labels = scales::number_format(accuracy = 0.1))  # Extend slightly above max
}

# histograms: 
hist_ba <- plot_hist(result_summary, metric = "balanced_accuracy", 
                     y_axis_label = "Count") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05))
  
hist_f1 <- plot_hist(result_summary, metric = "f1_score", 
                     y_axis_label = "Count",
                     x_axis_label = "F1 score") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0)) + 
  theme(axis.title.y = element_blank())

hist_f1_density <- plot_hist_density(result_summary, metric = "f1_score", 
                     y_axis_label = "Density",
                     x_axis_label = "F1 score") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0)) + 
  theme(axis.title.y = element_blank())


hist_f1a <- plot_hist(result_summary, metric = "f1_score", 
                     y_axis_label = "Count of instances",
                     x_axis_label = "F1 score") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0))

# Base‐R PDF device
# pdf(
#   file   = "hist_f1a.pdf",
#   width  = 5,    # inches
#   height = 4,
#   family = "Helvetica"   # or another installed font
# )
# print(hist_f1a)
# dev.off()     # close the file

combined_plot <- hist_ba + hist_f1 + 
  plot_layout(guides = "collect") +
  # Optionally, set the legend position (e.g., to the right or bottom)
  plot_annotation(theme = theme(legend.position = "right"))

combined_plot

hist_precision <- plot_hist(result_summary, metric = "precision", 
                            y_axis_label = "Count",
                            x_axis_label = "Precision") + 
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05)) +
  theme(axis.title.y = element_blank())

hist_recall <- plot_hist(result_summary, metric = "recall", 
                         y_axis_label = "Count",
                         x_axis_label = "Recall") + theme(axis.title.y = element_blank())

hist_specificity <- plot_hist(result_summary, metric = "specificity", 
                              y_axis_label = "Count",
                              x_axis_label = "Specificity") + theme(axis.title.y = element_blank())

hist_rmse <- plot_hist(result_summary, metric = "rmse", 
                       y_axis_label = "Count",
                       x_axis_label = "RMSE") + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1.0)) +
  theme(axis.title.y = element_blank())

hist_mse <- plot_hist(result_summary, metric = "mse", 
                      y_axis_label = "Count",
                      x_axis_label = "MSE") + 
  scale_x_continuous(labels = scales::number_format(accuracy = 1.0)) +
  theme(axis.title.y = element_blank())

combined_plot <- hist_precision + hist_recall + hist_specificity + hist_rmse + hist_mse +
  plot_layout(guides = "collect") +
  # Optionally, set the legend position (e.g., to the right or bottom)
  plot_annotation(theme = theme(legend.position = "right"))

combined_plot

# stats
# Run the t-test via formula interface
t_test_f1 <- t.test(f1_score ~ layer_comparison, 
                    data       = result_summary,
                    var.equal  = FALSE)  # Welch’s test

# 3. Print the full test
print(t_test_f1)
#> 
#> Welch Two Sample t-test
#> 
#> data:  f1_score by layer_comparison
#> t =  X.XXX, df =  Y.YY, p-value = Z.ZZZ
#> alternative hypothesis: true difference in means is not equal to 0
#> 95 percent confidence interval:
#>   LL UU
#> sample estimates:
#> mean in group Diagonal mean in group Off-diagonals 
#>                M₁                     M₂ 

# 4. Extract just the numbers you want
t_stat <- unname(t_test_f1$statistic)
df_val <- unname(t_test_f1$parameter)
p_val  <- t_test_f1$p.value

data.frame(
  t_value = t_stat,
  df      = df_val,
  p_value = p_val
)

# for poster:
custom_labels <- c("Diagonal" = "Single location",
                   "Off-diagonals" = "Location combination")

hist_f1p <- plot_hist(result_summary, metric = "f1_score", 
                      y_axis_label = "Count of instances",
                      x_axis_label = "F1 score") + 
  scale_fill_manual(values = custom_colors, labels = custom_labels) +
  labs(fill = "Prediction based on") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0)) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18))

# png(
#   filename = "hist_f1p.png",
#   width    = 8,           # width in inches
#   height   = 5,           # height in inches
#   units    = "in",        # could also be "px", "cm", etc.
#   res      = 300          # resolution in dots per inch
# )
# grid::grid.draw(hist_f1p)
# dev.off() 

# some stats

# Define the columns you want to summarize
eval_metrics <- c("specificity", "precision", "recall", "f1_score",
                  "balanced_accuracy", "mcc", "mse", "rmse")

# Summarize by layer_comparison
summary_stats <- result_summary %>%
  group_by(layer_comparison) %>%
  summarise(across(all_of(eval_metrics),
                   list(mean = mean,
                        sd   = sd,
                        max  = max,
                        min  = min),
                   .names = "{.col}_{.fn}"))

summary_stats

## ---- network size and density correlation with evaluators ----
# first we need to calculate the size and density of our networks
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict
# Load matrices. this is done differently for site scale and island scale.

### ---- island scale ----

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict
# Load matrices
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
A_l <- d$extended

# aggregate to island scale
# Extract numeric layer numbers
A_l <- A_l %>%
  mutate(layer_num = as.numeric(gsub("layer_", "", layer_from))) %>%
  mutate(aggregated_layer = ifelse(layer_num %% 2 == 1,
                                   paste0("layer_", layer_num, "_", layer_num + 1),
                                   paste0("layer_", layer_num - 1, "_", layer_num)))

# Aggregate data
aggregated_df <- A_l %>%
  group_by(aggregated_layer, node_from, node_to, type) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  mutate(layer_from = aggregated_layer, layer_to = aggregated_layer) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# Generate new layer names
unique_layers <- unique(aggregated_df$layer_from)  # Get unique aggregated layer names
new_layer_names <- paste0("layer_", seq_along(unique_layers))  # Generate new names (layer_1, layer_2, ...)

# Create a mapping table
layer_mapping <- data.frame(original_layer = unique_layers, new_layer = new_layer_names)

# Apply renaming in aggregated_df
aggregated_df <- aggregated_df %>%
  left_join(layer_mapping, by = c("layer_from" = "original_layer")) %>%
  mutate(layer_from = new_layer, layer_to = new_layer) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# View updated aggregated_df
print(aggregated_df)

# Total number of layers
num_layers <- length(unique(aggregated_df$layer_from))

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {

    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))

    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)

    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)

    node_to <- rownames(P) # for the results
    node_from <- colnames(P)

    # 1) find the species they have in common
    shared_plants      <- intersect(colnames(A), colnames(P))
    
    # 2) subset both matrices to exactly those shared species
    A <- A[ , shared_plants, drop = FALSE]
    P <- P[ , shared_plants, drop = FALSE]
    
    ### ---- creating a combined matrix C ----
    all_row_ids <- unique(c(rownames(A), rownames(P)))
    C <- matrix(0,
                nrow = length(all_row_ids),
                ncol = length(shared_plants),
                dimnames = list(all_row_ids, shared_plants))
    
    # Place A into C
    C[rownames(A), colnames(A)] <- A
    
    # Place P into C
    # Ensure that existing entries are not overwritten; sum overlapping entries
    C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                          NA, 
                                          C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    
    # Compute matrix properties
    nrow_A <- nrow(A)
    nrow_P <- nrow(P)
    nrow_C <- nrow(C)
    ncol_A <- ncol(A)
    ncol_P <- ncol(P)
    ncol_C <- ncol(C)
    size_A <- length(A)
    size_P <- length(P)
    size_C <- length(C)

    # Calculate density for A, P, and C
    density_A <- sum(A > 0) / length(A)
    density_P <- sum(P > 0) / length(P)
    density_C <- sum(C > 0) / length(C)

    # Add these values to the results table
    results <- rbind(results, data.frame(emln_id = emln_id,
                                         train_layer = layers_to_train,
                                         test_layer = layer_to_predict,
                                         nrow_A = nrow_A,
                                         ncol_A = ncol_A,
                                         size_A = size_A,
                                         density_A = density_A,   # Added density
                                         nrow_P = nrow_P,
                                         ncol_P = ncol_P,
                                         size_P = size_P,
                                         density_P = density_P,   # Added density
                                         nrow_C = nrow_C,
                                         ncol_C = ncol_C,
                                         size_C = size_C,
                                         density_C = density_C))  # Added density


  }
}

# View results
summary(results)

#results %>% write_csv('result_netsize_canaries_island_scale_50_itr.csv')

result_summary <- result_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

### ---- correlate netsize with evaluators ----

# if we want to use all of the results

df_long_1off <- result_summary %>%
  select(f1_score, balanced_accuracy, precision, recall, specificity, rmse, mse, size_P, density_P, size_C, density_C) %>%
  pivot_longer(
    cols = c(size_P, density_P, size_C, density_C),
    names_to = "measure_type",
    values_to = "measure_value"
  )

plot_netsize <- function(data, evaluator = "f1_score",
                         facet_labels = NULL,
                         evaluator_label = NULL,
                         title_text = NULL) {

  evaluator_sym <- rlang::sym(evaluator)  # Treat evaluator as a column

  # Calculate correlations
  cor_table <- data %>%
    group_by(measure_type) %>%
    summarise(
      cor_value = cor(!!evaluator_sym, measure_value, use = "complete.obs", method = "pearson"),
      p_value   = cor.test(!!evaluator_sym, measure_value, method = "pearson")$p.value,
      .groups = "drop"
    )

  # Create annotation table
  cor_table_annot <- cor_table %>%
    mutate(
      r_fmt = formatC(cor_value, format = "f", digits = 2),
      p_fmt = ifelse(
        p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),  # Scientific notation for very small p-values
        formatC(p_value, format = "f", digits = 3)   # Regular fixed format otherwise
      ),
      label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
    )


  # Build the plot
  plot <- ggplot(data, aes(x = measure_value, y = !!evaluator_sym)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "salmon") +
    facet_wrap(
      ~ measure_type,
      scales   = "free_x",
      labeller = as_labeller(facet_labels)
    ) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    geom_text(
      data    = cor_table_annot,
      aes(label = label_text),
      x       = Inf,
      y       = Inf,
      hjust   = 1.1,
      vjust   = 1.2,
      size    = 3.2,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Network feature",
      y = ifelse(is.null(evaluator_label), evaluator, evaluator_label),
      title = ifelse(is.null(title_text), paste(evaluator, "vs. network measures"), title_text)
    ) +
    theme_minimal() +
    tme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks = element_line(color = "black"),
      strip.text = element_text(size = 12)
    )

  return(plot)
}

netsize_island_f1 <- plot_netsize(
  data = df_long_1off,
  evaluator = "f1_score",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "F1 score"
)
netsize_island_f1

netsize_island_ba <- plot_netsize(
  data = df_long_1off,
  evaluator = "balanced_accuracy",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "Balanced accuracy"
)
netsize_island_ba

netsize_island_precision <- plot_netsize(
  data = df_long_1off,
  evaluator = "precision",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "Precision"
)
netsize_island_precision


netsize_island_recall <- plot_netsize(
  data = df_long_1off,
  evaluator = "recall",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "Recall"
)
netsize_island_recall

netsize_island_specificity <- plot_netsize(
  data = df_long_1off,
  evaluator = "specificity",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "Specificity"
)
netsize_island_specificity

netsize_island_rmse <- plot_netsize(
  data = df_long_1off,
  evaluator = "rmse",
  facet_labels = c(
    "size_C" = "Size of matrix C",
    "density_C" = "Density of matrix C",
    "size_P" = "Size of matrix P",
    "density_P" = "Density of matrix P"
  ),
  evaluator_label = "RMSE"
)
netsize_island_rmse

# netsize_island_mse <- plot_netsize(
#   data = df_long_1off,
#   evaluator = "mse",
#   facet_labels = c(
#     "size_C" = "Size of matrix C",
#     "density_C" = "Density of matrix C",
#     "size_P" = "Size of matrix P",
#     "density_P" = "Density of matrix P"
#   ),
#   evaluator_label = "MSE"
# )
# netsize_island_mse
# 
# 
plot_f1_rmse_vs_size_free_both <- function(data) {
  # build correlation table
  cor_table <- data %>%
    group_by(evaluator, measure_type) %>%
    summarise(
      cor_value = cor(evaluator_value, measure_value, use = "complete.obs"),
      p_value   = cor.test(evaluator_value, measure_value, method = "pearson")$p.value,
      .groups   = "drop"
    ) %>%
    mutate(
      r_fmt      = formatC(cor_value, format = "f", digits = 2),
      p_fmt      = ifelse(
        p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),
        formatC(p_value, format = "f", digits = 3)
      ),
      label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
    )

  ggplot(data, aes(x = measure_value, y = evaluator_value)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "salmon") +

    facet_grid(
      rows   = vars(evaluator),
      cols   = vars(measure_type),
      scales = "free",     # ← free both x and y per facet
      labeller = labeller(
        evaluator    = c(f1_score = "F1 score", rmse = "RMSE"),
        measure_type = c(size_P  = "Size of matrix P",
                         size_C  = "Size of matrix C")
      ),
      switch = "y"
    ) +

    geom_text(
      data        = cor_table,
      aes(label    = label_text),
      x           = Inf, y    = Inf,
      hjust       = 1.1, vjust = 1.2,
      size        = 3.2,
      inherit.aes = FALSE
    ) +

    scale_x_continuous(
      name   = "Network size",
      expand = expansion(mult = c(0.05, 0.1))
    ) +

    scale_y_continuous(
      name   = NULL,                # remove y title
      expand = expansion(mult = c(0.05, 0.1))
    ) +

    #labs(title = "F1 score and RMSE vs. Size of matrices P and C") +

    theme_minimal() +
    theme(
      strip.placement    = "outside",
      strip.text.x       = element_text(size = 14),
      strip.text.y.left  = element_text(size = 14, face = "bold", angle = 90),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks         = element_line(color = "black"),
      strip.background   = element_blank()
    )
}


df_f1_rmse_size <- result_summary %>%
  select(f1_score, rmse, size_P, size_C) %>%
  pivot_longer(cols = c(size_P, size_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f1_score, rmse), names_to = "evaluator", values_to = "evaluator_value")

netsize_f1_rmse <- plot_f1_rmse_vs_size_free_both(df_f1_rmse_size) + tme

# # pdf(
# #   file   = "netsize_f1_rmse.pdf",
# #   width  = 6,    # inches
# #   height = 6,
# #   family = "Helvetica"   # or another installed font
# # )
# # print(netsize_f1_rmse)
# # dev.off()     # close the file

## ---- degree impact and correlation with evaluators ----
### ---- calculate overall degree ----
# Step 1: Filter the data
df_filtered <- df %>%
  filter(itr == 1, original_links != 0)

# Step 2a: Calculate degree for each plant species (node_from)
plant_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(plant_degree = n(), .groups = "drop")

# Average plant degree by train_layer and test_layer
avg_plant_degree <- plant_degree %>%
  group_by(node_from) %>%
  summarise(avg_plant_degree = mean(plant_degree), .groups = "drop")

overall_plant_degree <- df_filtered %>% 
  group_by(node_from) %>% 
  summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")

overall_poll_degree <- df_filtered %>% 
  group_by(node_to) %>% 
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

# Step 2b: Calculate degree for each pollinator species (node_to)
pollinator_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(poll_degree = n(), .groups = "drop")

# Average pollinator degree by train_layer and test_layer
avg_pollinator_degree <- pollinator_degree %>%
  group_by(node_to) %>%
  summarise(avg_pollinator_degree = mean(poll_degree), .groups = "drop")

### ---- plot degree vs. number of never observed interactions ----
# here we examine if the algorithm assigns more links to species with higher degree.
# here by "island" we refer to a layer pair

df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

# Step 1: For each island and interaction, determine if the interaction was observed.
# Here we use `any(original_links == 1)` so that if the interaction is observed in at least one iteration, we count it.
df_island <- df %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )

table(df_island$observed)

# Step 2: Now, for each unique interaction, compute:
# - The proportion of islands where it was observed.
# - The average predicted probability (averaged over islands).
df_summary <- df_island %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

# 3. Remove interactions that were never observed and predicted as zero.
# df_summary <- df_summary %>% 
#   filter(!(avg_prop == 0 & avg_sigm_predicted < 0.5))

df_never_observed <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_never_observed_poll <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")

# for overall degree

never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")

never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")

# if we want to omit the Euphorbias
# never_plants_degree_overall <- never_plants_degree_overall %>%
#   filter(!node_from %in% c("Euphorbia_balsamifera_m", "Euphorbia_balsamifera_f"))

# correlation
df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 2)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, non-observed interactions",
    title = "Plants"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.2,   # Adjust depending on your data range
           label = label_text_plants,
           size = 5,
           color = "black")
plant_degree

df_to_correlate <- never_poll_degree_overall
df_to_correlate$x <- df_to_correlate$overall_poll_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_poll <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_poll
# Extract correlation coefficient and p-value
r_value <- round(correlation_poll$estimate, 2)
p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, /nnon-observed interactions",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.2,   # Adjust depending on your data range
           label = label_text_polls,
           size = 5,
           color = "black")
poll_degree

final_plot <- combine_plots(plant_degree, poll_degree)

# pdf(
#   file   = "degree.pdf",
#   width  = 8,
#   height = 5,
#   family = "Helvetica"
# )
# 
# grid::grid.draw(final_plot)
# 
# dev.off()

## ---- never-observed links ----
### ---- heatmap related to island proportion ----
# here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of islands in which these interactions were observed.

# order species by their degree
df_summary <- df_summary %>% left_join(overall_poll_degree, by="node_to")
df_summary <- df_summary %>% left_join(overall_plant_degree, by="node_from")

# Determine the order for plants based on overall_plant_degree
plant_order <- df_summary %>%
  distinct(node_from, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(node_from)

# Determine the order for pollinators based on overall_poll_degree
poll_order <- df_summary %>%
  distinct(node_to, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(node_to)

# Reset the levels for the species factors
df_summary$node_from <- factor(df_summary$node_from, levels = plant_order)
df_summary$node_to   <- factor(df_summary$node_to, levels = poll_order)


map_missing_links <- ggplot(df_summary, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop)) +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Proportion\nof islands\nobserved") +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
    aes(fill = avg_sigm_predicted),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2", 
                      name = "Average \npredicted \nprobability") +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_missing_links)

### ---- detect interactions that were never observed in the field yet consistently predicted to exist ----
# filter the interactions that were always observed as zeros yet predicted to exist
df_all_itr_zero <- df %>%
  group_by(node_from, node_to) %>%
  # Check if all rows for that interaction have original_links == 0
  filter(all(original_links == 0)) %>%
  filter(removed == 1) %>% 
  filter(predicted_prob_sigm > best_discrete_threshold) %>% 
  ungroup()

# filter only links that appear several times for the analysis
predicted_links <- df_all_itr_zero %>%
  group_by(node_from, node_to) %>%
  # Keep only those links/groups with >= 5 observations 
  filter(n() >= 5) %>%  
  summarise(
    n         = n(),
    mean_pred = mean(predicted_prob_sigm, na.rm = TRUE),
    sd_pred   = sd(predicted_prob_sigm, na.rm = TRUE),  # to check variance
    p_value   = if (sd_pred == 0) {
      NA_real_  # can't run a t-test if there's no variance
    } else {
      t.test(predicted_prob_sigm, mu = best_discrete_threshold)$p.value
    },
    .groups   = "drop"
  )

#write_csv(predicted_links, "predicted_non_observed_links_shared_species.csv")

most_probable_20 <- predicted_links %>%
  arrange(desc(mean_pred)) %>%
  head(20)

df_summary_sign <- predicted_links %>% filter(p_value < 0.05)

df_top_10 <- predicted_links %>%
  # Order by ascending p-value
  arrange(p_value) %>%
  # Take the first 10 rows
  slice(1:10)

#write.csv2(most_probable_20, "most_probable_20_links_shared_species.csv")

missing_links <- ggplot(df_top_10, 
       aes(x = reorder(paste(node_to, node_from, sep = " - "), mean_pred),
           y = mean_pred,
           fill = p_value < 0.05)) +
  geom_col(fill = "lightsteelblue") +
  geom_errorbar(aes(ymin = mean_pred - sd_pred, ymax = mean_pred + sd_pred),
                width = 0.2) +
  coord_flip() +
  #scale_fill_manual(name = "Significant?", values = c("gray70", "tomato")) +
  tme +
  scale_x_discrete(
    labels = function(x) sapply(x, function(lbl) {
      # 1) Replace underscores with a tilde (for spacing in plotmath)
      #    e.g., "Euphorbia_balsamifera_f" => "Euphorbia~balsamifera~f"
      lbl_tilde <- gsub("_", "~", lbl)
      # 2) Wrap in italic(), so the entire thing is in italics
      #    expression syntax: parse(text="italic(Euphorbia~balsamifera~f)")
      parse(text = paste0("italic(", lbl_tilde, ")"))
    })) +
  labs(
    x = "Link (pollinator - plant)",
    y = "Mean predicted value"
    #title = "Mean predicted value and significance"
  ) +
  theme(axis.title.y = element_text(margin = ggplot2::margin(r = 15)))

missing_links

### ---- missing links predicted by using shared species only vs. entire matrices ----
missing_links_shared_plants <- predicted_links
missing_links_all_species <- read.csv("predicted_non_observed_links.csv", row.names = NULL)

# Add an identifier to each row for comparison
missing_links_shared_plants$link_id <- paste(missing_links_shared_plants$node_from,
                                             missing_links_shared_plants$node_to,
                                              sep = "_")

missing_links_all_species$link_id <- paste(missing_links_all_species$node_from,
                                           missing_links_all_species$node_to,
                               sep = "_")

# Find shared and unique links
shared_links <- intersect(missing_links_shared_plants$link_id, missing_links_all_species$link_id)
unique_to_shared_species <- setdiff(missing_links_shared_plants$link_id, missing_links_all_species$link_id)
unique_to_missing_links <- setdiff(missing_links_all_species$link_id, missing_links_shared_plants$link_id)

# Percentage overlap
total_unique_links <- length(union(missing_links_all_species$link_id, missing_links_shared_plants$link_id))
overlap_percent <- length(shared_links) / total_unique_links * 100
print(paste0("Overlap: ", round(overlap_percent, 2), "%"))

# Combine the two tables and label their origin
shared_df <- rbind(
  missing_links_shared_plants %>% mutate(source = "shared_plants"),
  missing_links_all_species %>% mutate(source = "all_species")
)

# Classify each link into a category
shared_df <- shared_df %>%
  distinct(link_id, .keep_all = TRUE) %>%
  mutate(link_class = case_when(
    link_id %in% shared_links ~ "both_analyses",
    link_id %in% unique_to_shared_species ~ "shared_plants_only",
    link_id %in% unique_to_missing_links ~ "all_species_only"
  ))


ggplot(shared_df, aes(x = node_to, y = node_from)) +
  geom_tile(data = shared_df, aes(fill = link_class), alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "both_analyses" = "lightsteelblue",         # blue
      "shared_plants_only" = "plum3",    # purple
      "all_species_only" = "salmon1"      # orange
    ),
    name = "Link category",
    labels = c(
      "both_analyses" = "Both analyses",
      "shared_plants_only" = "Shared plants only",
      "all_species_only" = "All species only"
    )
  ) +
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  scale_y_discrete(
    labels = function(x) lapply(strsplit(x, "_"), function(y) {
      bquote(italic(.(paste(y, collapse = " "))))
    })
  ) + tme

# pie chart
link_summary <- tibble(
  category = c("Both analyses", "Shared plants only", "All species only"),
  count = c(length(shared_links),
            length(unique_to_shared_species),
            length(unique_to_missing_links))
) %>%
  mutate(percent = round(count / sum(count) * 100, 1),
         label = paste0(category, "\n", percent, "%"))

ggplot(link_summary, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    "Both analyses" = "lightsteelblue",        # blue
    "Shared plants only" = "thistle",   # purple
    "All species only" = "salmon1"     # orange
  )) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4,
            color = "white") +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("Link Overlap Between Tables") + tme

