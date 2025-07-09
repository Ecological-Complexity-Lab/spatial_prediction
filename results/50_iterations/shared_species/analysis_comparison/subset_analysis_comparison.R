# ---- comparing analyses quality ----
## ---- load libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(VennDiagram)
library(tidyverse)
library(pROC)
library(PRROC)
library(emln)
library(reshape2)
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
library(purrr)

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

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

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

## ---- load data ----
df_all <- read.csv('weighted__scaled_island_net_60_100_itr.csv')
df_shared_plants <- read.csv('weighted__scaled_island_net_60_100_itr_shared_plants.csv')
df_shared_pollinators <- read.csv('weighted_scaled_island_net_60_100_itr_shared_pollinators.csv')
df_shared_species <- read.csv('weighted__scaled_island_net_60_100_itr_shared_species.csv')

## ---- ok go ----
### ---- evaluation ----
analyze_predictions <- function(df, name = "Dataset") {
  thresholds <- seq(0, 1, by = 0.1)
  
  df_prepped <- df %>%
    filter(removed == 1) %>%
    mutate(
      predicted_prob   = sigmoid(predicted_values),
      original_binary  = if_else(original_links > 0, 1, 0)
    )
  
  df_thresh <- df_prepped %>%
    tidyr::expand_grid(threshold = thresholds) %>%
    mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
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
      rmse = sqrt(mse),
      .groups = "drop"
    ) %>%
    group_by(emln_id, train_layer, test_layer, threshold) %>%
    summarise(across(c(TP:rmse), mean, na.rm = TRUE), .groups = "drop")
  
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
  
  # Plot metric curves
  metric_plot <- ggplot(df_avg, aes(threshold, value, color = metric)) +
    geom_line(size = 1) +
    labs(title = paste(name, "- Average metrics"), x = "Threshold", y = "Value") +
    theme_minimal()
  
  df_wide <- df_avg %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    mutate(absdiff = abs(f1_score - balanced_accuracy)) %>%
    slice_min(absdiff, n = 1)
  
  best_threshold <- df_wide$threshold
  
  # Evaluate at optimal threshold
  df_removed <- df %>%
    filter(removed == 1) %>%
    mutate(
      predicted_prob_sigm = sigmoid(predicted_values),
      predicted_bin_sigm  = if_else(predicted_prob_sigm > best_threshold, 1, 0),
      original_binary      = if_else(original_links > 0, 1, 0)
    )
  
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
      nnse = 1 / (2 - nse),
      .groups = "drop"
    ) %>%
    group_by(emln_id, train_layer, test_layer) %>%
    summarise(across(TP:nnse, mean, na.rm = TRUE), .groups = "drop")
  
  # ROC and PR Curves
  auc_val <- plot_roc_curve(df_removed$original_binary, df_removed$predicted_values)
  plot_pr_curve(df_removed$original_binary, df_removed$predicted_values)
  
  list(
    best_threshold = best_threshold,
    result_summary = result_summary,
    metric_plot    = metric_plot,
    auc            = auc_val
  )
}


res_all <- analyze_predictions(df_all, "All")
res_shared_plants <- analyze_predictions(df_shared_plants, "Shared plants")
res_shared_pollinators <- analyze_predictions(df_shared_pollinators, "Shared pollinators")
res_shared_species <- analyze_predictions(df_shared_species, "Shared species")

### ---- roc curves ----
# 1. Put your results into a named list
roc_data_list <- list(
  All                 = df_all,
  Shared_plants       = df_shared_plants,
  Shared_pollinators  = df_shared_pollinators,
  Shared_species      = df_shared_species
)

# 2. For each subset, compute a pROC::roc object and turn it into a data.frame
roc_df <- purrr::imap_dfr(roc_data_list, function(df, subset_name) {
  df2 <- df %>%
    filter(removed == 1) %>%
    mutate(
      truth = as.numeric(original_links > 0),
      prob  = sigmoid(predicted_values)
    )
  roc_obj <- roc(df2$truth, df2$prob)
  
  # pull out the sensitivities / specificities
  tibble(
    subset     = subset_name,
    fpr        = 1 - roc_obj$specificities,
    tpr        =    roc_obj$sensitivities,
    threshold  =    roc_obj$thresholds,
    auc        = as.numeric(roc_obj$auc)  # same for every row
  )
})

# 3. Extract one AUC‐per‐subset for annotation
auc_labels <- roc_df %>%
  select(subset, auc) %>%
  distinct() %>%
  mutate(label = paste0("AUC = ", round(auc, 3)),
         x = 0.6,  # x/y are positions in FPR‐TPR space
         y = 0.2)

# 4. Plot with facets
# 4. plot, with custom colors + labeller
#    change these to whatever you like:
line_color   <- "lightsteelblue"
random_color <- "salmon"

ggplot(roc_df, aes(x = fpr, y = tpr)) +
  # your ROC line
  geom_line(color = line_color, size = 1) +
  # diagonal random‐guess line
  geom_abline(
    intercept = 0, slope = 1,
    linetype  = "dashed",
    color     = random_color
  ) +
  # facet and relabel
  facet_wrap(
    ~ subset, ncol = 2,
    labeller = labeller(subset = c(
      All                = "All",
      Shared_Plants      = "Shared plants",
      Shared_Pollinators = "Shared pollinators",
      Shared_Species     = "Shared species"
    ))
  ) +
  # annotate each panel with its AUC
  geom_text(
    data        = auc_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 0
  ) +
  coord_equal() +                     # square 0–1 axes
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    x = "False Positive Rate (1 − specificity)",
    y = "True Positive Rate (sensitivity)"
  ) +
  theme_minimal(base_size = 14) + tme

# pr curves

pr_data_list <- list(
  All                 = df_all,
  Shared_Plants       = df_shared_plants,
  Shared_Pollinators  = df_shared_pollinators,
  Shared_Species      = df_shared_species
)

pr_df <- imap_dfr(pr_data_list, function(df, subset_name) {
  df2 <- df %>%
    filter(removed == 1) %>%
    mutate(
      truth = as.numeric(original_links > 0),
      prob  = sigmoid(predicted_values)
    )
  
  pr_obj <- pr.curve(
    scores.class0 = df2$prob[df2$truth == 1],
    scores.class1 = df2$prob[df2$truth == 0],
    curve = TRUE
  )
  
  tibble(
    subset    = subset_name,
    recall    = pr_obj$curve[, 1],
    precision = pr_obj$curve[, 2],
    threshold = pr_obj$curve[, 3],
    aucpr     = pr_obj$auc.integral
  )
})

# 3. extract one AUPRC per subset for labeling
aucpr_labels <- pr_df %>%
  distinct(subset, aucpr) %>%
  mutate(
    label = paste0("AUPRC = ", round(aucpr, 3)),
    x = 0.2,  # adjust as needed
    y = 0.8
  )

# 4. choose your line color (reuse from ROC if you like)
line_color <- "lightsteelblue"

# 5. plot side-by-side PR curves
ggplot(pr_df, aes(x = recall, y = precision)) +
  # PR curve
  geom_line(color = line_color, size = 1) +
  # facet panels with pretty titles
  facet_wrap(
    ~ subset, ncol = 2,
    labeller = labeller(subset = c(
      All                = "All",
      Shared_Plants      = "Shared plants",
      Shared_Pollinators = "Shared pollinators",
      Shared_Species     = "Shared species"
    ))
  ) +
  # annotate each panel with its AUPRC
  geom_text(
    data        = aucpr_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 0
  ) +
  labs(
    x = "Recall",
    y = "Precision"
  ) +
  theme_minimal(base_size = 14) + tme

### ---- plot distances between analyses ----
df_all_combined <- bind_rows(
  res_all$result_summary             %>% mutate(dataset = "all"),
  res_shared_species$result_summary %>% mutate(dataset = "shared_species"),
  res_shared_plants$result_summary  %>% mutate(dataset = "shared_plants"),
  res_shared_pollinators$result_summary %>% mutate(dataset = "shared_pollinators")
)

# plot
# overall difference
# For nNSE
ggplot(df_all_combined, aes(x = dataset, y = nnse)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(method = "anova", label = "p.format",
                     label.y = max(df_all_combined$nnse, na.rm = TRUE) * 1.05) +
  labs(title = "nNSE across datasets", x = "Dataset", y = "nNSE") +
  theme_minimal() + tme

# For F1
ggplot(df_all_combined, aes(x = dataset, y = f1_score)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(method = "anova", label = "p.format",
                     label.y = max(df_all_combined$f1_score, na.rm = TRUE) * 1.05) +
  labs(x = "Dataset", y = "F1 Score") +
  theme_minimal() + tme

comparisons <- list(
  c("all","shared_species"), c("all","shared_plants"), c("all","shared_pollinators"),
  c("shared_species","shared_plants"), c("shared_species","shared_pollinators"), c("shared_plants","shared_pollinators")
)

# pairwise nnse
ggplot(df_all_combined, aes(x = dataset, y = nnse)) +
  geom_boxplot(notch = TRUE, alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.signif"
  ) +
  labs(x = "Network subset", y = "nNSE") +
  theme_minimal() + tme

# pairwise f1
ggplot(df_all_combined, aes(x = dataset, y = f1_score)) +
  geom_boxplot(notch = TRUE, alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.signif"
  ) +
  labs(x = "Network subset", y = "F1 score") +
  theme_minimal() + tme

### ---- missing interactions and venn ----

analyze_predicted_links <- function(df, best_threshold) {
  df %>%
    # Keep only interactions with all original_links == 0
    group_by(node_from, node_to) %>%
    filter(all(original_links == 0)) %>%
    ungroup() %>%
    
    # Keep only removed links
    filter(removed == 1) %>%
    
    # Add predicted probability if needed
    mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%
    
    # Keep only those above the threshold
    filter(predicted_prob_sigm > best_threshold) %>%
    
    # Keep only those appearing at least 5 times
    group_by(node_from, node_to) %>%
    filter(n() >= 5) %>%
    
    summarise(
      n         = n(),
      mean_pred = mean(predicted_prob_sigm, na.rm = TRUE),
      sd_pred   = sd(predicted_prob_sigm, na.rm = TRUE),
      p_value   = if (sd_pred == 0) NA_real_ else t.test(predicted_prob_sigm, mu = best_threshold)$p.value,
      .groups   = "drop"
    )
}

# Apply to each subset using its respective threshold
predicted_links_all <- analyze_predicted_links(df_all, res_all$best_threshold)
predicted_links_shared_species <- analyze_predicted_links(df_shared_species, res_shared_species$best_threshold)
predicted_links_shared_plants <- analyze_predicted_links(df_shared_plants, res_shared_plants$best_threshold)
predicted_links_shared_pollinators <- analyze_predicted_links(df_shared_pollinators, res_shared_pollinators$best_threshold)

# Combine 'node_to' and 'node_from' columns for each data frame to create interaction sets
all_species_interactions <- unique(paste(predicted_links_all$node_to, predicted_links_all$node_from, sep = "_"))
shared_species_interactions <- unique(paste(predicted_links_shared_species$node_to, predicted_links_shared_species$node_from, sep = "_"))
shared_plants_interactions <- unique(paste(predicted_links_shared_plants$node_to, predicted_links_shared_plants$node_from, sep = "_"))
shared_pollinators_interactions <- unique(paste(predicted_links_shared_pollinators$node_to, predicted_links_shared_pollinators$node_from, sep = "_"))

# Create a list of these sets to pass to the Venn diagram function
interaction_sets <- list(
  "All species" = all_species_interactions,
  "Shared species" = shared_species_interactions,
  "Shared plants" = shared_plants_interactions,
  "Shared pollinators" = shared_pollinators_interactions
)

# Plot the Venn diagram
venn.plot <- venn.diagram(
  x = interaction_sets,
  category.names = c("All species", "Shared species", "Shared plants", "Shared pollinators"),
  filename = NULL,  # You can save it as a file if needed
  output = TRUE,
  col = "transparent",  # Outline color of the circles
  fill = c("lightsteelblue", "orange", "lightseagreen", "thistle"),  # Colors for each circle
  alpha = 0.5,  # Transparency level for fill color
  cex = 1.5,  # Text size
  fontface = "bold",  # Text boldness
  fontfamily = "sans",  # Text font
  cat.cex = 1.5,  # Text size for the category names
  cat.fontface = "bold",  # Text font for category names
  cat.pos = 0,  # Positioning of the category names (0 means below)
  cat.dist = 0.1  # Distance of category names from circles
)

# Clear the plot window before drawing the Venn diagram
grid.newpage()  # Clear any existing plot

# Draw the Venn diagram
grid.draw(venn.plot)

### ---- netsizeee ----
#### ---- all species ----
res_all_summary <- res_all$result_summary

# first we need to calculate the size and density of our networks
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict

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

res_all_summary <- res_all_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

#### ---- shared species ----
res_species_summary <- res_shared_species$result_summary
results <- data.frame()

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
    shared_pollinators <- intersect(rownames(A), rownames(P))
    shared_plants      <- intersect(colnames(A), colnames(P))
    
    # 2) subset both matrices to exactly those shared species
    A <- A[shared_pollinators, shared_plants, drop = FALSE]
    P <- P[shared_pollinators, shared_plants, drop = FALSE]
    
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

res_species_summary <- res_species_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

#### ---- shared plants ----
res_plant_summary <- res_shared_plants$result_summary
results <- data.frame()
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

res_plant_summary <- res_plant_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

#### ---- shared pollinators ----
res_pollinator_summary <- res_shared_pollinators$result_summary
results <- data.frame()

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    # 1) find the pollinators they have in common
    shared_pollinators      <- intersect(rownames(A), rownames(P))
    
    # 2) subset both matrices to exactly those shared species
    A <- A[shared_pollinators, , drop = FALSE]
    P <- P[shared_pollinators, , drop = FALSE]
    
    ### ---- creating a combined matrix C ----
    all_col_ids <- unique(c(colnames(A), colnames(P)))
    C <- matrix(0,
                nrow = length(shared_pollinators),
                ncol = length(all_col_ids),
                dimnames = list(shared_pollinators, all_col_ids))
    
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

res_pollinator_summary <- res_pollinator_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

#### ---- plot netsize ----

# 1. Put your four summaries into a named list
summary_list <- list(
  All                = res_all_summary,
  `Shared plants`    = res_plant_summary,
  `Shared pollinators` = res_pollinator_summary,
  `Shared species`   = res_species_summary
)

combined <- imap_dfr(summary_list, ~ .x %>% mutate(subset = .y)) %>%
  select(subset, f1_score, rmse, nnse, size_P, size_C)

# 2) helper to pivot to long for any 'size_*' column
make_long <- function(df, size_col) {
  df %>%
    select(subset, all_of(size_col), f1_score, rmse, nnse) %>%
    pivot_longer(c(f1_score, rmse, nnse),
                 names_to  = "metric",
                 values_to = "value") %>%
    mutate(metric = recode(metric,
                           f1_score = "F1 score",
                           rmse     = "RMSE",
                           nnse     = "nNSE"
    ))
}

# 3) compute per‐facet correlation & p‐value + label coords
make_cor_labels <- function(df_long, size_col) {
  df_long %>%
    group_by(subset, metric) %>%
    summarise(
      cor_test = list(cor.test(.data[[size_col]], value)),
      .groups = "drop"
    ) %>%
    mutate(
      r     = map_dbl(cor_test, ~ .x$estimate),
      p     = map_dbl(cor_test, ~ .x$p.value),
      label = paste0("r = ",   round(r, 3),
                     "\np = ", signif(p, 2))
    ) %>%
    # pick a good spot for the label in each panel:
    left_join(
      df_long %>% group_by(subset, metric) %>%
        summarise(
          x = max(.data[[size_col]], na.rm=TRUE)*0.7,
          y = max(value, na.rm=TRUE)*0.9,
          .groups="drop"
        ),
      by = c("subset","metric")
    )
}

# ---- Plot 1: vs size_P ----
long_P     <- make_long(combined, "size_P")
labels_P   <- make_cor_labels(long_P, "size_P")

ggplot(long_P, aes(x = size_P, y = value)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +
  stat_cor(
    method       = "pearson",
    label.x.npc  = "right",
    label.y.npc  = "top",
    label.sep    = "\n",       # two‐line label: r on line1, p on line2
    size         = 3,
    hjust        = 1.05,       # nudge left just a touch
    vjust        = 1.05        # nudge down just a touch
  ) +
  facet_grid(metric ~ subset,
             scales = "free",    # free x & y per panel
             switch = "y") +
  labs(x = "Size of P", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 90),
    panel.grid.minor   = element_blank()
  ) + tme

# ---- netsize C ----
long_C     <- make_long(combined, "size_C")
labels_C   <- make_cor_labels(long_C, "size_C")

ggplot(long_C, aes(x = size_C, y = value)) +
  geom_point(color = "steelblue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "salmon") +
  stat_cor(
    method       = "pearson",
    label.x.npc  = "right",
    label.y.npc  = "top",
    label.sep    = "\n",       # two‐line label: r on line1, p on line2
    size         = 3,
    hjust        = 1.05,       # nudge left just a touch
    vjust        = 1.05        # nudge down just a touch
  ) +
  facet_grid(metric ~ subset,
             scales = "free",    # free x & y per panel
             switch = "y") +
  labs(x = "Size of C", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 90),
    panel.grid.minor   = element_blank()
  ) + tme

# ---- heatmap of per island interactions ----
df <- read.csv('weighted__scaled_island_net_60_100_itr.csv')

best_discrete_threshold <- res_all$best_threshold

df_filtered <- df %>%
  filter(itr == 1, original_links != 0)

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


overall_plant_degree <- df_filtered %>% 
  group_by(node_from) %>% 
  summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")

overall_poll_degree <- df_filtered %>% 
  group_by(node_to) %>% 
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

## Step 2: Now, for each unique interaction, compute:
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

#order species by their degree
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

## ---- for an island each time ----
# Function to generate a list of heatmaps, one per suffix level

plot_heatmaps_by_suffix <- function(df_island,
                                    overall_poll_degree,
                                    overall_plant_degree,
                                    threshold,
                                    suffix_levels = 1:7,
                                    tme             # your ggplot theme object
) {
  # 1) global summary to get ordering
  df_global <- df_island %>%
    group_by(node_from, node_to) %>%
    summarise(
      avg_prop           = mean(observed, na.rm = TRUE),
      avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
      n_islands          = n(),
      .groups = "drop"
    ) %>%
    left_join(overall_poll_degree,  by = "node_to") %>%
    left_join(overall_plant_degree, by = "node_from")
  
  plant_order <- df_global %>%
    distinct(node_from, overall_plant_degree) %>%
    arrange(desc(overall_plant_degree)) %>%
    pull(node_from)
  
  poll_order <- df_global %>%
    distinct(node_to, overall_poll_degree) %>%
    arrange(desc(overall_poll_degree)) %>%
    pull(node_to)
  
  plots <- vector("list", length(suffix_levels))
  names(plots) <- paste0("suffix_", suffix_levels)
  
  for (i in suffix_levels) {
    # subset & re-summarise
    df_sub_sum <- df_island %>%
      filter(str_detect(island_id, paste0("_", i, "$"))) %>%
      group_by(node_from, node_to) %>%
      summarise(
        avg_prop           = mean(observed, na.rm = TRUE),
        avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
        n_islands          = n(),
        .groups = "drop"
      ) %>%
      left_join(overall_poll_degree,  by = "node_to") %>%
      left_join(overall_plant_degree, by = "node_from") %>%
      # enforce the global factor levels here
      mutate(
        node_from = factor(node_from, levels = plant_order),
        node_to   = factor(node_to,   levels = poll_order)
      )
    
    # build the plot
    p <- ggplot(df_sub_sum, aes(x = node_to, y = node_from)) +
      # background layer
      geom_tile(aes(fill = avg_prop)) +
      scale_fill_gradient(
        low  = "white", high = "steelblue",
        name = "Proportion\nof islands\nobserved"
      ) +
      new_scale_fill() +
      # overlay layer
      geom_tile(
        data = df_sub_sum %>% filter(avg_prop == 0, avg_sigm_predicted > threshold),
        aes(fill = avg_sigm_predicted),
        alpha = 0.6
      ) +
      scale_fill_gradient(
        low  = "tan1", high = "tomato2",
        name = "Average\npredicted\nprobability"
      ) +
      # **force all species to appear** (even if no tiles for them)
      scale_x_discrete(limits = poll_order, drop = FALSE) +
      scale_y_discrete(
        limits = plant_order,
        labels = function(x) {
          lapply(strsplit(x, "_"), function(parts) {
            bquote(italic(.(paste(parts, collapse = " "))))
          })
        },
        drop = FALSE
      ) +
      theme_minimal() +
      labs(x = "Pollinator", y = "Plant") +
      theme(
        axis.text.x    = element_blank(),
        axis.text.y    = element_text(size = 8),
        legend.position = "bottom",
        legend.box      = "horizontal"
      ) +
      tme
    
    plots[[paste0("suffix_", i)]] <- p
  }
  
  plots
}

heatmap_list <- plot_heatmaps_by_suffix(
  df_island,
  overall_poll_degree,
  overall_plant_degree,
  threshold     = best_discrete_threshold,
  suffix_levels = 1:7,
  tme           = tme
)

# to print them:
for(p in heatmap_list) print(p)

