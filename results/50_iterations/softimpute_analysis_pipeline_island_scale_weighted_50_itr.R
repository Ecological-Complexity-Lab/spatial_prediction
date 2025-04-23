# ---- softImpute analysis pipeline for island scale ----
# this pipeline allows us to take the predictions of the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# for a first time run, run first site scale and then island scale to get the working dfs with evaluators. they are both needed for scale somparison.

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
  ggplot(data, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = FALSE, alpha = 0.4, color = "black") +
    theme_minimal() +
    labs(y = y_axis_label) +   # y-axis title is set via the function argument
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    tme + 
    scale_fill_manual(values = custom_colors) +
    stat_compare_means(method = "t.test", label = "p.signif", hide.ns = FALSE, 
                       comparisons = list(c("Diagonal", "Off-diagonals")),
                       label.y = stat_label_y,  # Adjust vertical position here
                       size = stat_size) +
    scale_y_continuous(limits = c (0.3, 0.9))
  
}

plot_hist <- function(data, metric, 
                      x_axis_label = "Balanced accuracy", 
                      y_axis_label = "Count") {
  ggplot(data, aes(x = .data[[metric]], fill = layer_comparison)) +
    geom_histogram(aes(y = ..count..), alpha = 0.4, color = "black", bins = 8, position = "dodge") +
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
                          left_label = "Number of predicted, non-observed links",
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
    rmse = sqrt(mse)
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
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup()

result_summary %>% write_csv('working_df_all_itr_60_weighted_scaled_island_50_itr.csv')
head(result_summary)
summary(result_summary)

### ---- distribution of evaluators ----
island_specificity <- ggplot(result_summary, aes(x = specificity)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Specificity",
       y = "Count") +
  tme

island_f1 <- ggplot(result_summary, aes(x = f1_score)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "F1 score",
       y = "Count") +
  tme

island_ba <- ggplot(result_summary, aes(x = balanced_accuracy)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Balanced accuracy",
       y = "Count") +
  tme

island_precision <- ggplot(result_summary, aes(x = precision)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Precision",
       y = "Count") +
  tme

island_recall <- ggplot(result_summary, aes(x = recall)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Recall",
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

combined_plots <- arrangeGrob(
  p1, p2, p3, p4, p5,
  ncol = 3, 
  nrow = 2
)
combined_with_axes <- arrangeGrob(
  combined_plots,
  #bottom = textGrob("F1 score", gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
  left   = textGrob("Count of instances", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

final_plot <- grid.arrange(
  combined_with_axes,
  ncol = 3,
  widths = c(2, 0.3, 0.3)
)

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
                        y_axis_label = "F1 score", stat_label_y = 0.85, stat_size = 6)

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

# histograms: 
hist_ba <- plot_hist(result_summary, metric = "balanced_accuracy", 
                     y_axis_label = "Count")
hist_f1 <- plot_hist(result_summary, metric = "f1_score", 
                     y_axis_label = "Count",
                     x_axis_label = "F1 score")

hist_f1 <- hist_f1 + theme(axis.title.y = element_blank())

combined_plot <- hist_ba + hist_f1 + 
  plot_layout(guides = "collect") +
  # Optionally, set the legend position (e.g., to the right or bottom)
  plot_annotation(theme = theme(legend.position = "right"))

combined_plot
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
    
    ### ---- creating a combined matrix C ----
    all_row_ids <- unique(c(rownames(A), rownames(P)))
    all_col_ids <- unique(c(colnames(A), colnames(P)))
    C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                dimnames = list(all_row_ids, all_col_ids))
    
    # Place A into C
    C[rownames(A), colnames(A)] <- A
    
    # Place P into C
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

results %>% write_csv('result_netsize_canaries_island_scale_50_itr.csv')

result_summary <- result_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

### ---- correlate netsize with evaluators ----

# if we want to use all of the results
df_long_1off <- result_summary %>%
  pivot_longer(
    cols = c(size_P, density_P, size_C, density_C),
    names_to = "measure_type",
    values_to = "measure_value"
  )

# For each variable, compute correlation with f1_score:
cor_table <- df_long_1off %>%
  group_by(measure_type) %>%
  summarise(
    cor_value = cor(f1_score, measure_value, use = "complete.obs", method = "pearson"),
    p_value   = cor.test(f1_score, measure_value, method = "pearson")$p.value
  ) %>%
  ungroup()

cor_table

cor_table_annot <- cor_table %>%
  mutate(
    # round correlation to 3 decimals, no scientific notation
    r_fmt  = formatC(cor_value, format = "f", digits = 2),
    # round p-value to 4 decimals, no scientific notation
    p_fmt  = formatC(p_value,  format = "f", digits = 4),
    label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
  )

# Create a named vector for renaming facets
facet_labels <- c(
  "size_C" = "Size of matrix C",
  "density_C" = "Density of matrix C",
  "size_P" = "Size of matrix P",
  "density_P" = "Density of matrix P"
)

netsize_island <- ggplot(df_long_1off, aes(x = measure_value, y = f1_score)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +
  facet_wrap(
    ~ measure_type,
    scales   = "free_x",
    labeller = as_labeller(facet_labels)  # Use the named vector
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
    color   = "black"
  ) +
  labs(
    x = "Network feature",
    y = "F1 score",
    title = "F1 vs. network measures - all data points"
  ) +
  theme_minimal() +
  tme +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black")
  )
print(netsize_island)

## ---- Jaccard correlation with evaluators ----

# Initialize a data frame to store combined results for all layer combinations
results_jaccard <- data.frame()

### ---- island scale ----
# Total number of layers
num_layers <- length(unique(aggregated_df$layer_from)) # we have it from netsize calculation

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
    P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
    
    # 1) Jaccard pollinators
    poll_train <- rownames(A)[ rowSums(A) > 0 ]
    poll_test  <- rownames(P)[ rowSums(P) > 0 ]
    intersection_poll <- length(intersect(poll_train, poll_test))
    union_poll        <- length(union(poll_train, poll_test))
    jaccard_poll <- if (union_poll == 0) NA else intersection_poll / union_poll
    
    # 2) Jaccard plants
    plants_train <- colnames(A)[ colSums(A) > 0 ]
    plants_test  <- colnames(P)[ colSums(P) > 0 ]
    intersection_plants <- length(intersect(plants_train, plants_test))
    union_plants        <- length(union(plants_train, plants_test))
    jaccard_plants <- if (union_plants == 0) NA else intersection_plants / union_plants
    
    # 3) Jaccard edges
    pairs_train <- which(A > 0, arr.ind = TRUE)
    pairs_train_strings <- apply(pairs_train, 1, function(rc) {
      paste(rownames(A)[rc[1]], colnames(A)[rc[2]], sep = "_")
    })
    
    pairs_test <- which(P > 0, arr.ind = TRUE)
    pairs_test_strings <- apply(pairs_test, 1, function(rc) {
      paste(rownames(P)[rc[1]], colnames(P)[rc[2]], sep = "_")
    })
    intersection_edges <- length(intersect(pairs_train_strings, pairs_test_strings))
    union_edges        <- length(union(pairs_train_strings, pairs_test_strings))
    jaccard_edges <- if (union_edges == 0) NA else intersection_edges / union_edges
    
    # store the result
    results_jaccard <- rbind(
      results_jaccard,
      data.frame(
        emln_id = emln_id,
        train_layer = layers_to_train,
        test_layer  = layer_to_predict,
        jaccard_pollinators = jaccard_poll,
        jaccard_plants      = jaccard_plants,
        jaccard_edges       = jaccard_edges
      )
    )
  }
}


head(results_jaccard)
result_summary <- result_summary %>%
  left_join(results_jaccard, by = c("train_layer", "test_layer")) # add to results table

#### ---- plot ----
make_facet_scatter_plot <- function(data,
                                    evaluator = "f1_score", 
                                    pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                    names_to = "jaccard_type", 
                                    values_to = "jaccard_value",
                                    x_lab = "Jaccard similarity",
                                    y_lab = "F1 score",
                                    plot_title = "F1 vs. Jaccard - All data (site)",
                                    facet_scales = "free_x") {
  
  # Reshape data from wide to long format for the specified pivot columns
  df_long <- data %>% 
    pivot_longer(cols = all_of(pivot_cols), 
                 names_to = names_to, 
                 values_to = values_to)
  
  # For each facet (jaccard_type), compute correlation between the evaluator and jaccard_value
  cor_table <- df_long %>%
    group_by(!!sym(names_to)) %>%
    summarise(
      cor_value = cor(.data[[evaluator]], .data[[values_to]], use = "complete.obs", method = "pearson"),
      p_value   = cor.test(.data[[evaluator]], .data[[values_to]], method = "pearson")$p.value
    ) %>%
    ungroup()
  
  # Create annotations with formatted correlation coefficients and p-values
  cor_table_annot <- cor_table %>%
    mutate(
      r_fmt      = formatC(cor_value, format = "f", digits = 2),
      p_fmt      = formatC(p_value, format = "f", digits = 3),
      label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
    )
  
  # Construct the faceted scatter plot
  plot <- ggplot(df_long, aes_string(x = values_to, y = evaluator)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "thistle") +
    facet_wrap(as.formula(paste("~", names_to)), scales = facet_scales) +
    scale_x_continuous(labels = number_format(accuracy = 0.1)) +
    # Place the correlation annotation in the upper-right corner of each facet
    geom_text(data = cor_table_annot,
              aes(label = label_text),
              x = Inf,
              y = Inf,
              hjust = 1.1,
              vjust = 1.2,
              size = 3.2,
              color = "black") +
    labs(x = x_lab, y = y_lab, title = plot_title) +
    theme_minimal() +
    tme +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(color = "black")
    )
  
  return(plot)
}

# all data points
all_island <- make_facet_scatter_plot(data = result_summary, 
                                       evaluator = "f1_score",
                                       pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                       x_lab = "Jaccard similarity",
                                       y_lab = "F1 score",
                                       plot_title = "F1 vs. Jaccard - all data (island)",
                                       facet_scales = "free_x")

all_island

# off-diagonals
canary_results_diags_isl <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer != test_layer)

offs_island <- make_facet_scatter_plot(data = canary_results_diags_isl, 
                                     evaluator = "f1_score",
                                     pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                     x_lab = "Jaccard similarity",
                                     y_lab = "F1 score",
                                     plot_title = "F1 vs. Jaccard - off-diagonals (island)",
                                     facet_scales = "free_x")

offs_island

## ---- partner fidelity correlation with evaluators ----

# filter out cases in which train = test layer
df_fidelity <- df %>% filter(train_layer == test_layer) %>% 
  filter(original_links != 0) %>% filter(itr == 1)

### ---- plants fidelity ----
# 1. For each plant (node_from) and layer, gather the pollinators (node_to).
#    (Assuming 'train_layer' is the relevant layer ID—adapt as needed.)
df_plant_partners <- df_fidelity %>%
  distinct(node_from, train_layer, node_to) %>%
  group_by(node_from, train_layer) %>%
  summarise(partners = list(unique(node_to)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen <- df_plant_partners %>%
  group_by(node_from) %>%
  summarise(
    mean_sorensen_plants = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

df_merged <- df %>%
  left_join(df_sorensen, by = "node_from")

### ---- pollinators fidelity ----

df_pollinators <- df_fidelity %>%
  distinct(node_to, train_layer, node_from) %>%
  group_by(node_to, train_layer) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen_pollinators <- df_pollinators %>%
  group_by(node_to) %>%
  summarise(
    mean_sorensen_pollinators = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

df_merged <- df_merged %>%
  left_join(df_sorensen_pollinators, by = "node_to")

### --- correlation with evaluators ----
summary_df <- df_merged %>%
  group_by(train_layer, test_layer) %>%
  summarise(
    avg_sorensen_plants = mean(mean_sorensen_plants, na.rm = TRUE),
    avg_sorensen_pollinators = mean(mean_sorensen_pollinators, na.rm = TRUE)
  ) %>%
  ungroup()

result_summary <- result_summary %>% 
  left_join(summary_df, by = c("train_layer", "test_layer"))
  
working_df_offs <- result_summary %>% filter (train_layer != test_layer)

correlation <- cor.test(working_df_offs$f1_score, working_df_offs$avg_sorensen_pollinators, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

pollinator_fidelity_cor <- ggplot(working_df_offs, aes(x = avg_sorensen_pollinators, y = f1_score)) +
  geom_point(color = "thistle", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "F1 score",
       title = "Island scale (pollinators)") +
  tme +
  annotate("text",
              x = Inf,
              y = Inf,
              hjust = 1.1,
              vjust = 1.2,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black")+
  theme(axis.title.y = element_blank()) # for the unified plot

correlation <- cor.test(working_df_offs$f1_score, working_df_offs$avg_sorensen_plants, use = "complete.obs", method = "pearson")
correlation

# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

plant_fidelity_cor <- ggplot(working_df_offs, aes(x = avg_sorensen_plants, y = f1_score)) +
  geom_point(color = "darkseagreen3", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "F1 score",
       title = "Island scale (plants)") +
  tme +
  annotate("text", x = Inf, y = Inf, label = label_text,
            hjust = 1.1, vjust = 1.1, size = 3.5, color = "black") +
  theme(axis.title.y = element_blank()) # for the unified plot


grid.arrange(
  arrangeGrob(plant_fidelity_cor, pollinator_fidelity_cor, ncol = 2),
  left = textGrob("F1 score", rot = 90, gp = gpar(fontsize = 13, fontface = "bold"))
)
  
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

# correlation
df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 3)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, non-observed interactions",
    title = "Plants: island scale"
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
r_value <- round(correlation_poll$estimate, 3)
p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, non-observed interactions",
    title = "Pollinators: island scale"
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


ggplot(df_summary, aes(x = node_to, y = node_from)) +
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

predicted_links %>%
  arrange(desc(mean_pred)) %>%
  head(10)

#write.csv(predicted_links, 'links_probability_to_be_non_zeros.csv')

df_summary_sign <- predicted_links %>% filter(p_value < 0.05)

df_top_10 <- predicted_links %>%
  # Order by ascending p-value
  arrange(p_value) %>%
  # Take the first 10 rows
  slice(1:10)

ggplot(df_top_10, 
       aes(x = reorder(paste(node_to, node_from, sep = " - "), mean_pred),
           y = mean_pred,
           fill = p_value < 0.05)) +
  geom_col(fill = "steelblue") +
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
  )

## ---- distance effect ----
### ---- add distances and location names ----
distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL)

# proceed for both island scale and site scale and compare the trends
### ---- island scale ----
# Function to extract island names (removes "_site_X")
extract_island <- function(name) {
  gsub("_site_[12]", "", name)
}

# Create new table with averaged distances at the island level
distance_island_table <- distance_table %>%
  mutate(
    from_island = extract_island(from),
    to_island = extract_island(to)
  ) %>%
  group_by(from_island, to_island) %>%
  summarise(
    avg_distance_m = mean(distance_m),
    avg_distance_km = mean(distance_km),
    .groups = "drop"
  ) %>%
  mutate(
    avg_distance_m = ifelse(from_island == to_island, 0, avg_distance_m),
    avg_distance_km = ifelse(from_island == to_island, 0, avg_distance_km)
  ) %>%
  rename(from = from_island, to = to_island)  # Rename after calculation

# Print result
print(distance_island_table)

# Modify the 'from' and 'to' columns in distance_island_table
distance_island_table <- distance_island_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

# Add names and distances to the main table
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Step 1: Create a new grouped tibble
new_layer_names <- net_name %>%
  mutate(group_id = (layer_id + 1) %/% 2) %>%  # Group pairs into 1, 2, 3...
  group_by(group_id) %>%
  summarise(name = gsub(" site.*", "", first(name)), .groups = "drop")  # Keep only location name

# Add to main table
result_summary <- result_summary %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)

result_summary_island <- result_summary

# Add to main table
result_summary_island <- result_summary_island %>%
  left_join(
    distance_island_table,
    by = c("train_layer_name" = "from", "test_layer_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_layer_name == test_layer_name,
                               0,              # distance = 0 if same site
                               avg_distance_km))   # otherwise, keep joined distance

### ---- site scale ----
# Add names and distances to the main table
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Modify the 'from' and 'to' columns in distance_table
distance_table <- distance_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

# Assuming your lookup tibble is called net_name and has columns layer_id and name
#result_summary_site <- read_csv('working_df_site_scaled_evaluators_distance_50_itr_net_60.csv')
result_site <- read_csv('canary_weighted_scaled_site_net_60_50_itr.csv')

df_removed_site <- result_site %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

result_summary_site <- df_removed_site %>%
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
    rmse = sqrt(mse)
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
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup()


result_summary_site <- result_summary_site %>%
  # Join to add train_layer_name
  left_join(net_name %>% 
              rename(train_layer = layer_id, 
                     train_layer_name = name), 
            by = "train_layer") %>%
  # Join to add test_layer_name
  left_join(net_name %>% 
              rename(test_layer = layer_id, 
                     test_layer_name = name), 
            by = "test_layer")


# Add to main table
result_summary_site <- result_summary_site %>%
  left_join(
    distance_table,
    by = c("train_layer_name" = "from", "test_layer_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_layer_name == test_layer_name,
                               0,              # distance = 0 if same site
                               distance_km))   # otherwise, keep joined distance

## ---- distance correlation with evaluators ----
make_cor_plot <- function(data, evaluator, 
                          distance_col = "distance_km", 
                          x_lab = "Geographical distance (km)",
                          y_lab = NULL,
                          extra_theme = NULL) {
  # Use evaluator as y_lab if no alternative is provided
  if (is.null(y_lab)) {
    y_lab <- evaluator
  }
  
  # Compute correlation between evaluator and distance
  correlation <- cor.test(data[[evaluator]], data[[distance_col]], 
                          use = "complete.obs", method = "pearson")
  r_value <- round(correlation$estimate, 3)
  p_value <- formatC(correlation$p.value, format = "f", digits = 4)
  label_text <- paste0("r = ", r_value, ", p = ", p_value)
  
  # Create plot with label in the upper right corner using Inf coordinates
  plot <- ggplot(data, aes_string(x = distance_col, y = evaluator)) +
    geom_point(color = "salmon2", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +
    labs(x = x_lab, y = y_lab) +
    # The following places the label at the upper right of the plot area
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.1, size = 3.5, color = "black")
  
  # Optionally add additional theme modifications
  if (!is.null(extra_theme)) {
    plot <- plot + extra_theme
  }
  
  return(plot)
}

cor_plot_site <- make_cor_plot(result_summary_site, evaluator = "f1_score", extra_theme = tme)
cor_plot_isl  <- make_cor_plot(result_summary_island, evaluator = "f1_score", extra_theme = tme)

# To combine the plots:
p1 <- cor_plot_site + 
  ggtitle("Site scale") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 0.3), "cm"))
p2 <- cor_plot_isl +
  ggtitle("Island scale") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 0.3), "cm"))

combined_plots <- arrangeGrob(
  p1, p2,
  ncol = 2,
  widths = c(1, 1)
)
combined_with_axes <- arrangeGrob(
  combined_plots,
  bottom = textGrob("Geographical distance (km)", 
                    gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
  left   = textGrob("F1 score", rot = 90, 
                    gp = gpar(fontsize = 14, fontface = "bold"))
)
final_plot <- grid.arrange(
  combined_with_axes,
  ncol = 2,
  widths = c(2, 0.3)
)

final_plot

## ---- plot heatmaps ----
island_heatmap_recall <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = recall)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "Recall") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_recall)

island_heatmap_f1 <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "F1 score") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_f1)

island_heatmap_ba <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "Balanced \naccuracy") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_ba)

island_heatmap_specificity <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = specificity)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "Specificity") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_specificity)

island_heatmap_precision <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = precision)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "Precision") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(island_heatmap_precision)

# # For an evaluator "recall" on a data frame 'result_summary_island'
# # (assuming tme is a predefined ggplot theme you want to apply)
# island_heatmap_ba <- make_heatmap_plot(data = result_summary_island, 
#                                            eval_col = "balanced_accuracy", 
#                                            train_col = "train_layer_name", 
#                                            test_col = "test_layer_name",
#                                            x_lab = "Added layer",
#                                            y_lab = "Predicted layer",
#                                            fill_lab = "Balanced \naccuracy",
#                                            extra_theme = tme)
# 
# island_heatmap_f1 <- make_heatmap_plot(data = result_summary_island, 
#                                        eval_col = "f1_score", 
#                                        train_col = "train_layer_name", 
#                                        test_col = "test_layer_name",
#                                        x_lab = "Added layer",
#                                        y_lab = "Predicted layer",
#                                        fill_lab = "F1 score",
#                                        extra_theme = tme)
# 
# island_heatmap_recall <- make_heatmap_plot(data = result_summary_island, 
#                                            eval_col = "recall", 
#                                            train_col = "train_layer_name", 
#                                            test_col = "test_layer_name",
#                                            x_lab = "Added layer",
#                                            y_lab = "Predicted layer",
#                                            fill_lab = "Recall",
#                                            extra_theme = tme)
# 
# island_heatmap_precision <- make_heatmap_plot(data = result_summary_island, 
#                                            eval_col = "precision", 
#                                            train_col = "train_layer_name", 
#                                            test_col = "test_layer_name",
#                                            x_lab = "Added layer",
#                                            y_lab = "Predicted layer",
#                                            fill_lab = "Precision",
#                                            extra_theme = tme)
# 
# island_heatmap_specificity <- make_heatmap_plot(data = result_summary_island, 
#                                               eval_col = "specificity", 
#                                               train_col = "train_layer_name", 
#                                               test_col = "test_layer_name",
#                                               x_lab = "Added layer",
#                                               y_lab = "Predicted layer",
#                                               fill_lab = "Specificity",
#                                               extra_theme = tme)
# 
# p1 <- island_heatmap_f1 +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))
# 
# p2 <- island_heatmap_ba +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))
# 
# p3 <- island_heatmap_recall +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))
# 
# p4 <- island_heatmap_precision +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))
# 
# p5 <- island_heatmap_specificity +
#   theme(legend.position = "none",
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))
# 
# combined_plots <- arrangeGrob(
#   p1, p2, p3, p4, p5,
#   ncol = 3, 
#   nrow = 2
# )
# combined_with_axes <- arrangeGrob(
#   combined_plots,
#   #bottom = textGrob("F1 score", gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
#   left   = textGrob("Predicted location", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
# )
# 
# final_plot <- grid.arrange(
#   combined_with_axes,
#   ncol = 3,
#   widths = c(2, 0.3, 0.3)
# )
## ---- compare scales ----
# make a long list
df1_labeled <- result_summary_site %>%
  mutate(scale = "Site")

df2_labeled <- result_summary_island %>%
  mutate(scale = "Island")

df_combined <- bind_rows(df1_labeled, df2_labeled)

df_long <- df_combined %>%
  pivot_longer(
    cols = c("f1_score", "recall", "precision", "balanced_accuracy", "mcc", "specificity"),
    names_to = "metric",
    values_to = "value"
  )

metrics <- c("f1_score", "recall", "precision", "balanced_accuracy", "mcc", "specificity")

results <- lapply(metrics, function(metric) {
  test_normality_site <- shapiro.test(result_summary_site[[metric]])$p.value
  test_normality_island <- shapiro.test(result_summary_island[[metric]])$p.value
  
  if (test_normality_site > 0.05 & test_normality_island > 0.05) {
    test <- t.test(result_summary_site[[metric]], result_summary_island[[metric]], var.equal = FALSE)
  } else {
    test <- wilcox.test(result_summary_site[[metric]], result_summary_island[[metric]])
  }
  
  data.frame(
    Metric = metric,
    Test = ifelse(test_normality_site > 0.05 & test_normality_island > 0.05, "T-test", "Wilcoxon"),
    P_value = test$p.value
  )
})

results_df <- do.call(rbind, results)
print(results_df)

# Define significance function
get_pvalue_asterisks <- function(p) {
  if (p < 0.001) return("***")  # Highly significant
  else if (p < 0.01) return("**")  # Very significant
  else if (p < 0.05) return("*")  # Significant
  else return("ns")  # Not significant
}

stat_results <- lapply(metrics, function(metric) {
  data_metric <- df_long %>% filter(metric == !!metric)  # Filter for the specific metric
  
  test <- t.test(value ~ scale, data = data_metric)  # Perform t-test
  
  p_value <- test$p.value
  significance <- get_pvalue_asterisks(p_value)
  
  data.frame(
    metric = metric,
    p_value = p_value,
    significance = significance
  )
})

stat_results_df <- do.call(rbind, stat_results)

# Merge significance levels with the dataset
df_long <- df_long %>%
  left_join(stat_results_df, by = "metric")

# Create the boxplot with significance annotations
ggplot(df_long, aes(x = metric, y = value, fill = scale)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Comparison of Performance", x = "Metric", y = "Value") +
  scale_fill_manual(values = c("Site" = "lightsteelblue2", "Island" = "wheat2")) +  # Custom colors
  scale_x_discrete(labels = c(
    "f1_score" = "F1 score",
    "recall" = "Recall",
    "precision" = "Precision",
    "balanced_accuracy" = "Balanced \naccuracy",
    "mcc" = "MCC",
    "specificity" = "Specificity"
  )) +  # Properly formatted labels
  stat_compare_means(aes(group = scale), method = "t.test", label = "p.signif", 
                     label.y = max(df_long$value, na.rm = TRUE) + 0.05,
                     size = 5)  + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14), ) + tme

## ---- variable importance ----
# analyze only one off-diagonal
df_off <- result_summary_island %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer)
# subset relavant columns
df_subset <- df_off %>% select(f1_score, distance_km,	avg_sorensen_plants,	avg_sorensen_pollinators,	jaccard_pollinators,	jaccard_plants,	jaccard_edges, size_P,	density_P, size_C,	density_C)

### ---- autocorrelation check ----
# 2. Compute correlation matrix
cor_mat <- cor(df_subset, use = "complete.obs")
#write.csv(cor_mat, "island_autocorrelation.csv")

# 3. Visualize (optional)

corrplot(cor_mat, 
         method = "number",      # or "circle", "color", etc.
         type = "upper",         # upper/lower/full
         tl.cex = 0.7,           # text label size
         number.cex = 0.7,       # correlation coefficient size
         tl.col = "black"       # text label color (optional)
)

### ---- linear regression ----
# Fit a linear model predicting f1_score from all other numeric predictors
lm_fit <- lm(f1_score ~ distance_km +	avg_sorensen_plants +	avg_sorensen_pollinators +	jaccard_pollinators +	jaccard_plants +	jaccard_edges + size_P +	density_P + size_C +	density_C,
             data = df_subset)

summary(lm_fit)

# Random Forest (only with numeric columns)
rf_fit <- randomForest(f1_score ~ ., data = df_subset, importance = TRUE)

# Check variable importance
importance(rf_fit)
varImpPlot(rf_fit)

# Extract importance
imp <- importance(rf_fit) 
# For regression: imp is a matrix with columns: %IncMSE, IncNodePurity
# For classification: imp often has two columns per measure.
# We'll assume %IncMSE and IncNodePurity are present.

# Turn it into a data frame for easier plotting
imp_df <- as.data.frame(imp)
imp_df$f1_score <- rownames(imp_df)  # Keep variable names in a column

# Example for %IncMSE
ggplot(imp_df, aes(x = reorder(f1_score, `%IncMSE`), y = `%IncMSE`)) +
  geom_bar(stat = "identity", fill = "lightsteelblue") +
  coord_flip() +
  labs(x = "F1 score", 
       y = "% Increase in MSE") +
  theme_minimal() + tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

## ---- pca of latent traits ----
