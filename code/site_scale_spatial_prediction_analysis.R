# ---- Site scale: predicting interactions across space with SVD ----
# this pipeline allows us to predict missing links at the site level using the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# here we focus on site scale, which is shown in Appendix S1.
# stages are according to the pipeline figure (Fig. 1).
### code for publication ###
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
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
library(stringr)
library(softImpute)
library(ecodist)
library(rstatix)
library(ggrepel)

source("code/common.R") 

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- parameters ----
prop_ones_to_remove <- 0.2 # proportion of existing links to withhold
n_sim <- 50 # number of random link withholding and prediction iterations
set.seed(42) # the answer to everything

## ---- functions ----
# building matrices for combining matrices, calculating network size and density
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

# predict links using softImpute
implement_impute <- function(C, k, lambda) {
  # Apply softImpute
  
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  
  # Debias the fit to remove regularization effects
  # fit <- deBias(C, fit)
  
  # Reconstruct the matrix
  C_reconstructed <- softImpute::complete(C, fit)
  
  # back-transform C to original scale (was centered using biScale)
  C_reconstructed_orig <- C_reconstructed +
    outer(row_centers[rownames(C)], col_centers[colnames(C)], "+")
  
  # Extract the reconstructed P matrix from C_reconstructed
  P_reconstructed <- C_reconstructed_orig[rownames(P), colnames(P)]
  
  # Combine indices of removed ones and zeros
  if (is.null(dim(remove_indices))) { # handles when remove_indices has only one row
    test_indices <- rbind(
      data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  } else {
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = rep(1, nrow(remove_indices))),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
    )
  }
  
  # Get the row and column names of the test links
  test_rows <- rownames(P)[test_indices$row]  # These are the "node_to"
  test_cols <- colnames(P)[test_indices$col]  # These are the "node_from"
  
  # Actual labels and predictions
  original_links <- P_original[cbind(test_rows, test_cols)]
  predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
  
  # Store the results with node information
  results <- data.frame(k = k,
                        lambda = lambda,
                        original_links = original_links,
                        predicted_values = predicted_values,
                        node_to = test_rows,
                        node_from = test_cols,
                        removed = 1 # mark these links as removed
  )
  
  # ---- (A) Enumerate ALL edges in P
  all_edges <- expand.grid(
    node_to = rownames(P),
    node_from = colnames(P),
    k = k,
    lambda = lambda,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Fill in the original link values from P_original
  all_edges$original_links <- mapply(
    function(r, c) P_original[r, c],
    all_edges$node_to,
    all_edges$node_from
  )
  
  # Helper data frame of removed edges
  removed_edges_idx <- data.frame(
    node_to = test_rows,
    node_from = test_cols,
    stringsAsFactors = FALSE
  )
  
  # ---- (B) Subset edges NOT removed
  not_removed <- all_edges[
    !paste(all_edges$node_to, all_edges$node_from) %in%
      paste(removed_edges_idx$node_to, removed_edges_idx$node_from), # all node pairs that are not in th removed list
  ]
  not_removed$removed <- 0
  not_removed$predicted_values <- NA
  
  # ---- (C) Combine removed + not removed
  not_removed$k <- k
  not_removed$lambda <- lambda
  
  list_results <- list(results=results, not_removed=not_removed)
  
  return(list_results)
}

# plotting functions
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

make_facet_scatter_plot <- function(data,
                                    evaluator = "f05_score", 
                                    pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                    names_to = "jaccard_type", 
                                    values_to = "jaccard_value",
                                    x_lab = "Jaccard similarity",
                                    y_lab = "F0.5 score",
                                    plot_title = NULL,
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
      r_fmt = formatC(cor_value, format = "f", digits = 2),
      p_fmt = ifelse(
        p_value < 0.001,
        formatC(p_value, format = "e", digits = 2),  # Scientific notation
        formatC(p_value, format = "f", digits = 3)   # Otherwise
      ),
      label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)
    )
  
  # Facet labels (renaming)
  facet_labels <- c(
    jaccard_edges = "Interaction overlap",
    jaccard_plants = "Plants overlap",
    jaccard_pollinators = "Pollinators overlap"
  )
  
  panel_labels <- tibble::tibble(
    !!names_to := pivot_cols,
    panel_label = c("(c)", "(b)", "(a)")
  )
  
  # Construct the faceted scatter plot
  plot <- ggplot(df_long, aes_string(x = values_to, y = evaluator)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "thistle") +
    facet_wrap(as.formula(paste("~", names_to)), 
               scales = facet_scales,
               labeller = as_labeller(facet_labels)) +
    scale_x_continuous(labels = number_format(accuracy = 0.02)) +
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
    geom_text(
      data = panel_labels,
      aes(label = panel_label),
      x = -Inf,
      y = Inf,
      hjust = -0.4,
      vjust = 1.4,
      size = 5,
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    theme_minimal() +
    tme +
    theme(
      strip.text = element_text(size = 12),  # <-- Facet titles larger and bold
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(color = "black"),
      axis.text.x     = element_text(size = 12),   
      axis.text.y     = element_text(size = 12))
  
  return(plot)
}

make_simple_correlation_plot <- function(data,
                                         x_var,
                                         evaluator,
                                         x_lab = NULL,
                                         y_lab = NULL,
                                         plot_title,
                                         point_color,
                                         trend_color = "steelblue",
                                         x_axis_blank = FALSE,
                                         y_axis_blank = FALSE) {
  
  # Perform correlation test
  correlation <- cor.test(data[[evaluator]], data[[x_var]], use = "complete.obs", method = "pearson")
  
  # Extract correlation coefficient and p-value
  r_value <- round(correlation$estimate, 2)
  p_value <- ifelse(
    correlation$p.value < 0.001,
    formatC(correlation$p.value, format = "e", digits = 2),  # scientific for very small
    formatC(correlation$p.value, format = "f", digits = 3)   # fixed format otherwise
  )  
  
  label_text <- paste0("r = ", r_value, ", p = ", p_value)
  
  # Build the plot
  p <- ggplot(data, aes_string(x = x_var, y = evaluator)) +
    geom_point(color = point_color, alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = trend_color) +
    labs(
      x = x_lab,
      y = y_lab,
      title = plot_title
    ) +
    tme +
    annotate("text",
             x = Inf, y = Inf,
             hjust = 1.1, vjust = 1.2,
             label = label_text,
             size = 4,
             color = "black")
  
  # Optionally remove axis titles
  if (y_axis_blank) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank())
  }
  if (x_axis_blank) {
    p <- p + theme(axis.title.x = element_blank())
  }
  
  return(p)
}

# Final master function to make the full double plot
make_full_correlation_plot <- function(data,
                                       evaluator = "f05_score",
                                       pollinator_x = "avg_sorensen_pollinators",
                                       plant_x = "avg_sorensen_plants",
                                       shared_x_lab = "Mean Sorensen similarity",
                                       shared_y_lab = NULL) {
  
  if (is.null(shared_y_lab)) {
    # If user doesn't specify left y-axis label, use the evaluator name nicely formatted
    shared_y_lab <- gsub("_", " ", evaluator)
    shared_y_lab <- str_to_title(shared_y_lab)
  }
  
  # Create plots (with no x labels)
  pollinator_plot <- make_simple_correlation_plot(
    data = data,
    x_var = pollinator_x,
    evaluator = evaluator,
    x_lab = NULL,  # No individual x-label
    y_lab = NULL,  # No individual y-label
    plot_title = "Pollinators",
    point_color = "thistle",
    trend_color = "steelblue",
    x_axis_blank = TRUE,
    y_axis_blank = TRUE
  )
  
  plant_plot <- make_simple_correlation_plot(
    data = data,
    x_var = plant_x,
    evaluator = evaluator,
    x_lab = NULL,  # No individual x-label
    y_lab = NULL,  # No individual y-label
    plot_title = "Plants",
    point_color = "darkseagreen3",
    trend_color = "steelblue",
    x_axis_blank = TRUE,
    y_axis_blank = FALSE
  )
  
  plant_plot <- plant_plot + theme(plot.margin = ggplot2::margin(4, 10, 4, 10))
  pollinator_plot <- pollinator_plot + theme(plot.margin = ggplot2::margin(4, 10, 4, 10))
  
  
  # Arrange plots side by side
  # plots_side_by_side <- arrangeGrob(
  #   plant_plot, pollinator_plot,
  #   ncol = 2
  # )
  # Arrange plots side by side with equal widths
  plots_side_by_side <- arrangeGrob(
    plant_plot, pollinator_plot,
    ncol = 2,
    widths = unit.c(unit(1.13, "null"), unit(1, "null"))  # Equal widths
  )
  
  
  # Add shared axis labels
  final_plot <- grid.arrange(
    plots_side_by_side,
    left = textGrob(shared_y_lab, rot = 90, gp = gpar(fontsize = 13, fontface = "bold")),
    bottom = textGrob(shared_x_lab, gp = gpar(fontsize = 13, fontface = "bold"))
  )
  
  return(final_plot)
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

make_cor_plot <- function(data, evaluator, 
                          distance_col = "distance_km", 
                          x_lab = "Geographic distance (km)",
                          y_lab = NULL,
                          extra_theme = NULL) {
  # Use evaluator as y_lab if no alternative is provided
  if (is.null(y_lab)) {
    y_lab <- evaluator
  }
  
  # Compute correlation between evaluator and distance
  correlation <- cor.test(data[[evaluator]], data[[distance_col]], 
                          use = "complete.obs", method = "pearson")
  r_value <- round(correlation$estimate, 2)
  p_value <- ifelse(
    correlation$p.value < 0.001,
    formatC(correlation$p.value, format = "e", digits = 2),  # scientific for very small
    formatC(correlation$p.value, format = "f", digits = 3)   # fixed format otherwise
  ) 
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


## ---- 1. prediction ----
# Load matrices
d <- load_emln(emln_id)
A_l <- d$extended

# Extract numeric layer numbers
A_l <- A_l %>%
  mutate(layer_num = as.numeric(gsub("layer_", "", layer_from))) %>%
  mutate(aggregated_layer = ifelse(layer_num %% 2 == 1, 
                                   paste0("layer_", layer_num, "_", layer_num + 1),
                                   paste0("layer_", layer_num - 1, "_", layer_num)))

# Total number of layers
num_layers <- length(unique(A_l$layer_from))

# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

# Loop through all combinations of layers_to_train and layer_to_predict
for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    ### ---- a. withhold links in P ----
    # map out the 0s and 1s in P
    num_1_to_remove <- floor(sum(P>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
    ones_in_P <- which(P > 0, arr.ind = TRUE)
    
    num_0_to_remove <- num_1_to_remove
    prop_0_removed <- num_0_to_remove / sum(P == 0, na.rm = T)
    zeros_in_P <- which(P == 0, arr.ind = TRUE)
    
    # debug print
    print(paste("1 remove:", num_1_to_remove))
    print(paste("all 1   :", nrow(ones_in_P)))
    print(paste("0s to remove:", num_0_to_remove))
    print(paste("all zeros   :", nrow(zeros_in_P)))
    print(paste("prop of zeros removed   : ", prop_0_removed))
    
    # Randomly select zeros to withhold - bootstrapping
    bootstrapping_results <- NULL
    P_original <- P # save it for later
    
    for (i in 1:n_sim) {
      # remove 1s
      remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
      P[remove_indices] <- NA  # Set removed links to NA
      
      # sample 0s
      zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
      P[zeros_to_remove_indices] <- NA
      
      ### ---- creating a combined matrix C ----
      # Combine A and P into a single matrix C with NAs representing missing data
      all_row_ids <- unique(c(rownames(A), rownames(P)))
      all_col_ids <- unique(c(colnames(A), colnames(P)))
      C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                  dimnames = list(all_row_ids, all_col_ids))
      
      # Place A into C
      C[rownames(A), colnames(A)] <- A
      
      # Place P into C
      if (layers_to_train != layer_to_predict){
        # Ensure that existing entries are not overwritten; sum overlapping entries
        C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                              NA, 
                                              C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
      } else {
        # if this predicts using the same layer, don't sum it to itself
        C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                              NA, 
                                              (C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])/2)
      }
      

      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
      # save centers before overwriting (for back-transforming later)
      row_centers <- attr(C, "biScale:row")$center      # named vector, length = nrow(C)
      col_centers <- attr(C, "biScale:column")$center   # named vector, length = ncol(C)
      sum(is.na(C))
      
      ### ---- b. prediction with SVD ----
      k_values <- c(2)
      lam0 <- lambda0(C)
      lambda_values <- c(lam0)
      
      # Initialize variables to store the best results
      results <- data.frame(k = integer(),
                            lambda = numeric(),
                            original_links = numeric(),
                            predicted_values = numeric())
      not_removed_all <- NULL
      
      # Loop over all combinations of k and lambda
      for (k in k_values) {
        for (lambda in lambda_values) {
          # imputation
          r <- implement_impute(C, k, lambda)
          
          results <- rbind(results, r$results)
          not_removed_all <- rbind(not_removed_all, r$not_removed)
        }
      }
      
      ### ---- save results for current k/lambda combination ----
      # After finishing the k/lambda loops, append the 'results' to 'combined_results'
      # ---- (D) Append to combined_results
      complete_edges_all <- rbind(results, not_removed_all)
      complete_edges_all$itr <- i
      bootstrapping_results <- rbind(bootstrapping_results, complete_edges_all)
      
      # reset P
      P <- P_original
    }
    
    combined_results <- rbind(
      combined_results,
      cbind(
        data.frame(
          emln_id = emln_id,
          train_layer = layers_to_train,
          test_layer = layer_to_predict,
          prop_ones_removed = prop_ones_to_remove,
          amount_of_removed_1 = num_1_to_remove,
          amount_of_removed_0 = num_0_to_remove,
          prop_0_removed = prop_0_removed
        ),
        bootstrapping_results
      )
    )
  }
}
# combined_results includes predictions for all combinations of islands, 50 iterations of links withholding and prediction for each combination

## ---- 2. analysis ----
summary(combined_results)

# convert negatives to zeros
df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# know thy network - what species do we have in the system?
plant_species <- unique(df$node_from)           # get unique names
pollinator_species <- unique(df$node_to)

### ---- c. evaluation ----
### ---- selecting optimal threshold ----
# select the threshold for classifying links as 1s or 0s based on max f0.5

# 0) set an array of thresholds
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
  tidyr::expand_grid(threshold = thresholds) %>%  
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity = TN / (TN + FP),
    precision   = TP / (TP + FP),
    recall      = TP / (TP + FN),
    f05_score   = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
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
    f05_score = mean(f05_score, na.rm = TRUE),
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
      f05_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold,
               names_to  = "metric",
               values_to = "value")


# 1) pivot to wide so F0.5 and balanced_accuracy are columns
df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# 2) find the threshold with the optimal f0.5 score
best_discrete <- df_wide %>%
  slice_max(f05_score, n = 1)

# results:
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold

df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

### ---- create evaluation table ----
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
    f05_score = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
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
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()

head(result_summary)
summary(result_summary) # result_summary includes evaluation results across all iterations for each combination of islands 

### ---- Fig. S3: distribution of evaluators with/without external data ----
# this analysis shows us if predictions made using added information from other locations (off-diagonals in layer-to-layer predictions, as a heatmap) is any better than not adding any information (cases on the diagonal)
result_summary <- result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Single location",
    train_layer != test_layer ~ "Added location"
  ))

custom_colors <- c("Single location" = "steelblue",
                   "Added location" = "thistle")

# plot the histogram: 

hist_f05a <- plot_hist(result_summary, metric = "f05_score", 
                     y_axis_label = "Count of instances",
                     x_axis_label = "F0.5 score") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.05))

hist_f05a

# # Base‐R PDF device
pdf(
  file   = "results/paper_figs/hist_f05_site.pdf",
  width  = 5,    # inches
  height = 4,
  family = "Helvetica"   # or another installed font
)
print(hist_f05a)
dev.off()     # close the file

# stats
# run t-test via formula interface
t_test_f05 <- t.test(f05_score ~ layer_comparison, 
                    data       = result_summary,
                    var.equal  = FALSE)  # Welch’s test

# 3. Print the full test
print(t_test_f05)

# do we need welch/wilcoxon?
# first normality check
result_summary %>%
  group_by(layer_comparison) %>%
  shapiro_test(f05_score)

# variance check
result_summary %>% levene_test(f05_score ~ layer_comparison)
# all is good, we can use t-test.

# 4. Extract just the numbers you want
t_stat <- unname(t_test_f05$statistic)
df_val <- unname(t_test_f05$parameter)
p_val  <- t_test_f05$p.value

data.frame(
  t_value = t_stat,
  df      = df_val,
  p_value = p_val
)

### ---- e. ecological inference ----
### ---- network size and density correlation with evaluators ----
# here we calculate the size and density of the networks and correlate them with evaluation metrics.

#### ---- calculate the size and density of our networks ----
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

# Loop through all combinations of emln_id, layers_to_train, and layer_to_predict
# Load matrices
d <- load_emln(emln_id)
A_l <- d$extended

# Total number of layers
num_layers <- length(unique(A_l$layer_from))

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
    node_to <- rownames(P) # for the results
    node_from <- colnames(P)
    
    ### creating a combined matrix C 
    all_row_ids <- unique(c(rownames(A), rownames(P)))
    all_col_ids <- unique(c(colnames(A), colnames(P)))
    C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                dimnames = list(all_row_ids, all_col_ids))
    
    # Place A into C
    C[rownames(A), colnames(A)] <- A

    # Place P into C
    if (layers_to_train != layer_to_predict){
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    } else {
      # if this predicts using the same layer, don't sum it to itself
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            (C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])/2)
    }
    
    
    # make them all binary for count
    A[A>0] <- 1
    P[P>0] <- 1
    C[C>0] <- 1
    
    # Compute matrix properties
    nrow_A <- nrow(A)
    nrow_P <- nrow(P)
    nrow_C <- nrow(C)
    ncol_A <- ncol(A)
    ncol_P <- ncol(P)
    ncol_C <- ncol(C)
    size_A <- nrow(A) + ncol(A)
    size_P <- nrow(P) + ncol(P)
    size_C <- nrow(C) + ncol(C)
    
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

# view results
summary(results)


# add them to the results table
result_summary <- result_summary %>%
  left_join(results, by = c("train_layer", "test_layer")) # add to results table

# summerize
# first add layer names
# add distances to the main table
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

result_summary <- result_summary %>%
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

# now we can summarise
df_summary <- result_summary %>%
  group_by(test_layer_name) %>%
  summarise(
    size_P   = mean(size_P, na.rm = TRUE),
    density_P  = mean(density_P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(size_P))  # change to asc() for smallest first

df_summary

overall_sd_density <- sd(df_summary$density_P, na.rm = TRUE)

#### ---- Fig. S6: correlate network size with evaluators ----

df_netsize <- result_summary %>%
  select(f05_score, nnse, size_P, density_P, size_C, density_C) %>%
  pivot_longer(
    cols = c(size_P, density_P, size_C, density_C),
    names_to = "measure_type",
    values_to = "measure_value"
  )


df_f05_nnse_size <- result_summary %>%
  select(f05_score, nnse, size_P, size_C) %>%
  pivot_longer(cols = c(size_P, size_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f05_score, nnse), names_to = "evaluator", values_to = "evaluator_value")

netsize_f05_nnse <- plot_f05_nnse_vs_size_free_both(df_f05_nnse_size) + tme # Fig. 5
netsize_f05_nnse

pdf(
  file   = "results/paper_figs/site_netsize_f05_nnse.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(netsize_f05_nnse)
dev.off()     # close the file

#### ---- Fig. S7: density ----

df_f05_nnse_density <- result_summary %>%
  select(f05_score, nnse, density_P, density_C) %>%
  pivot_longer(cols = c(density_P, density_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f05_score, nnse), names_to = "evaluator", values_to = "evaluator_value")

netdensity_f05_nnse <- plot_f05_nnse_vs_density_free_both(df_f05_nnse_density) + tme
netdensity_f05_nnse

pdf(
  file   = "results/paper_figs/site_netdensity_f05_nnse.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(netdensity_f05_nnse)
dev.off()     # close the file

### ---- Fig. S2: heatmap ----
lims <- range(result_summary$f05_score, na.rm = TRUE)
mid_val <- mean(lims)

site_heatmap_f05 <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f05_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "salmon2", 
                       midpoint = mid_val, na.value = "gray") +  # Set NA values to gray
  labs(x = "Auxiliary site", y = "Target site", fill = "F0.5 score") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(site_heatmap_f05)

pdf(
  file   = "results/paper_figs/site_heatmap_f05.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(site_heatmap_f05)
dev.off()     # close the file

### ---- Fig. S4: Jaccard correlation with evaluators ----
#### ---- calculate Jaccard ----
# here we calculate the Jaccard index for every pair of networks
# Initialize a data frame to store combined results for all layer combinations
results_jaccard <- data.frame()

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)
    P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
    
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

#### ---- plot Jaccard ----
# we filter only pairs of different islands for this analysis, since Jaccard index for the same island is 1 for all islands
canary_results_jaccard <- result_summary %>%
  filter(train_layer != test_layer)

jaccard_site_f05 <- make_facet_scatter_plot(data = canary_results_jaccard, 
                                        evaluator = "f05_score",
                                        pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                        x_lab = "Jaccard similarity",
                                        y_lab = "F0.5 score",
                                        facet_scales = "free_x")
jaccard_site_f05

# # Base‐R PDF device
pdf(
  file   = "results/paper_figs/jaccard_site_f05.pdf",
  width  = 7,    # inches
  height = 3.5,
  family = "Helvetica"   # or another installed font
)
print(jaccard_site_f05)
dev.off()     # close the file

## ---- distance decay ----
#### ---- add distances and location names ----
distance_table <- read.csv("data/distance_between_sites_canary.csv", row.names = NULL)

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

result_summary_site <- result_summary %>%
  select(-train_layer_name, -test_layer_name) %>%   # remove old columns if they exist
  left_join(
    net_name %>%
      rename(train_layer = layer_id, train_layer_name = name),
    by = "train_layer"
  ) %>%
  left_join(
    net_name %>%
      rename(test_layer = layer_id, test_layer_name = name),
    by = "test_layer"
  )

# Add to main table
result_summary_site <- result_summary_site %>%
  left_join(
    distance_table,
    by = c("train_layer_name" = "from", "test_layer_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_layer_name == test_layer_name,
                               0,              # distance = 0 if same site
                               distance_km))   # otherwise, keep joined distance


# remove sites form within the same island - only use information from different islands for distance decay
result_summary_site_dif <- result_summary_site %>% filter(train_layer != test_layer)
# extract island names
result_summary_site_dif$train_island <- sub("^(\\w+).*", "\\1", result_summary_site_dif$train_layer_name)
result_summary_site_dif$test_island <- sub("^(\\w+).*", "\\1", result_summary_site_dif$test_layer_name)

# filter rows where island names are different
filtered_results <- result_summary_site_dif[result_summary_site_dif$train_island != result_summary_site_dif$test_island, ]

### ---- Fig. S5: plot distance decay ----
cor_plot_site_dif_f05 <- make_cor_plot(filtered_results, evaluator = "f05_score", extra_theme = tme) + labs(y = "F0.5 score")
pdf(
  file   = "results/paper_figs/cor_plot_site_dif_f05.pdf",
  width  = 4,
  height = 4,
  family = "Helvetica"
)
print(cor_plot_site_dif_f05)
dev.off()     # close the file

#### ---- MRM for site scale ----
layers_site <- sort(unique(c(result_summary_site_dif$train_layer, result_summary_site_dif$test_layer)))

# initialize empty matrices
f05_mat_site      <- matrix(NA, nrow=length(layers_site), ncol=length(layers_site),
                           dimnames=list(layers_site, layers_site))
dist_mat_km_site <- f05_mat_site

# fill in each cell [i,j] with the corresponding f05_score and distance_km
for(i in layers_site) for(j in layers_site) {
  # subset rows where train=i and test=j
  sub <- result_summary_site_dif[result_summary_site_dif$train_layer==i & result_summary_site_dif$test_layer==j, ]
  if(nrow(sub)==1) {
    f05_mat_site[i,j]      <- sub$f05_score
    dist_mat_km_site[i,j] <- sub$distance_km
  }
}

# because MRM uses symmetric distance matrices, average [i,j] & [j,i]

f05_sym_site      <- sym_average(f05_mat_site)
dist_sym_km_site <- sym_average(dist_mat_km_site)

# convert to “dist” objects (lower triangle)
dist_f05_site      <- as.dist(f05_sym_site)
dist_km_site     <- as.dist(dist_sym_km_site)

# run the MRM
# — this will regress the F0.5‐distance matrix on the geographic–distance matrix
set.seed(42)   # for reproducibility of permutations
mrm_out_site <- MRM(dist_f05_site ~ dist_km_site, nperm=999)

# results
print(mrm_out_site)

### ---- Fig. S1: compare scales ----
# read island scale results
combined_results_island <- 
  readRDS(file = paste0("results/predictions_island_scale.rds"))
combined_results_island <- combined_results_island %>% 
  filter(k == 2) %>% 
  filter(!(input_lambda %in%  c(1, 5, 50, 100)))

df_island <- combined_results_island %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# create the evaluation table
df_removed_island <- df_island %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

result_summary_island <- df_removed_island %>%
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
    f05_score = mean(f05_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()


df_long <- bind_rows(
  result_summary_site   %>% mutate(scale = "Site"),
  result_summary_island %>% mutate(scale = "Island")
) %>%
  pivot_longer(
    cols      = c("f05_score", "nnse"),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c(
    "f05_score", "nnse"
  )))

# test for assumptions
# normality

df_long %>%
  group_by(metric, scale) %>%
  shapiro_test(value)   # Shapiro-Wilk normality test

# for island scale, data is not normally distributed (F0.5 scores), so better use wilcoxon 
# plot
# pretty facet titles with units:
metric_labels <- c(
  f05_score         = "F0.5 score",
  nnse              = "NNSE"
)

panel_labels <- tibble::tibble(
  metric = names(metric_labels),
  panel_label = c("(a)", "(b)")
)

nnse_f05_scales <- ggplot(df_long, aes(x = scale, y = value, fill = scale)) +
  geom_boxplot(
    notch        = TRUE,
    outlier.size = 1,
    position     = position_dodge(width = 0.75)
  ) +
  facet_wrap(
    ~ metric,
    scales        = "free_y",
    labeller      = as_labeller(metric_labels),
    ncol          = 4,
    switch        = "y"
  ) +
  stat_compare_means(
    method         = "wilcox.test",
    label          = "p.format",
    p.format.args  = list(
      digits     = 2,
      scientific = TRUE
    ),
    label.y        = Inf,
    vjust          = 1.5,
    label.x        = 1.45,
    tip.length     = 0.01,
    size           = 3.5
  ) +
  geom_text(
    data = panel_labels,
    aes(label = panel_label),
    x = -Inf,
    y = Inf,
    hjust = -0.4,
    vjust = 1.4,
    size = 5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Site"   = "lightsteelblue2",
                               "Island" = "wheat2")) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement       = "outside",
    strip.text.y.left     = element_text(
      angle = 90,
      face  = "bold",
      size  = 12
    ),
    axis.text.x           = element_blank(),
    axis.ticks.x          = element_blank(),
    legend.position       = "bottom"
  ) +
  tme

nnse_f05_scales

# #Base‐R PDF device
pdf(
  file   = "results/paper_figs/nnse_f05_scales.pdf",
  width  = 7,    # inches
  height = 4,
  family = "Helvetica"   # or another installed font
)
print(nnse_f05_scales)
dev.off()     # close the file










