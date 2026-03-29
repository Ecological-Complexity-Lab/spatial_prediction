# ---- Predicting interactions across space with SVD ----
# this pipeline allows us to predict missing links using the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# here we focus on island scale, but there is a code for site-scale analysis and comparison between scales.
# stages are according to the pipeline figure (Fig. 1).
### code for publication ###
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)

# this is for installing the EMLN package (Frydman et al. 2023): designed for handling and analysing of ecological multilayer networks;
# in this pipeline it is used to import published multilayer network data
package.list=c("tidyverse", "magrittr","igraph","Matrix","DT","hablar","devtools")
loaded <-  package.list %in% .packages()
package.list <-  package.list[!loaded]
installed <-  package.list %in% .packages(TRUE)
if (!all(installed)) install.packages(package.list[!installed],repos="http://cran.rstudio.com/")

# Install EMLN only if not already installed
if (!requireNamespace("emln", quietly = TRUE)) {
  devtools::install_github("Ecological-Complexity-Lab/emln")
}

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
library(pROC)
library(PRROC)
library(magick)
library(pdftools)

source("code/common.R")


## ---- parameters ----
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

# function to extract island names (removes "_site_X")
extract_island <- function(name) {
  gsub("_site_[12]", "", name)
}


# plotting

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
      text = element_text(size = 18),
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
    theme_minimal() +
    tme +
    theme(
      strip.text = element_text(size = 12),  # <-- Facet titles larger and bold
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(color = "black"),
      axis.text.x     = element_text(size = 9)                          )
  
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
  
  return(combined_with_axes)
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

combine_two_plots <- function(p1, p2,
                              x_axis_label = "Geographic distance (km)",
                              y_axis_label = "F0.5 score",
                              p1_title = "Site scale",
                              p2_title = "Island scale",
                              margins = unit(c(0.5, 0.5, 1, 0.3), "cm"),
                              axis_title_fontsize = 14,
                              axis_title_fontface = "bold") {
  
  # Adjust first plot
  p1 <- p1 +
    ggtitle(p1_title) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      plot.margin = margins
    )
  
  # Adjust second plot
  p2 <- p2 +
    ggtitle(p2_title) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      plot.margin = margins
    )
  
  # Combine p1 and p2 side by side
  combined_plots <- arrangeGrob(
    p1, p2,
    ncol = 2,
    widths = c(1.1, 1)
  )
  
  # Add global x and y axis labels
  combined_with_axes <- arrangeGrob(
    combined_plots,
    bottom = textGrob(
      x_axis_label,
      gp = gpar(fontsize = axis_title_fontsize, fontface = axis_title_fontface),
      vjust = -1.5
    ),
    left = textGrob(
      y_axis_label,
      rot = 90,
      gp = gpar(fontsize = axis_title_fontsize, fontface = axis_title_fontface)
    )
  )
  
  # Final arrangement
  final_plot <- grid.arrange(
    combined_with_axes,
    ncol = 2,
    widths = c(2, 0.01)
  )
  
  return(final_plot)
}

## ---- 1. prediction ----
### ---- load matrices ----
# load and mold data
aggregated_df <- load_and_mold_data_for_prediction(emln_id)

# save aggregated network to a file
# write.csv(aggregated_df, file = "prediction_pipeline_for_publication/results/network_island_scale.csv", row.names = FALSE)

# Total number of layers
num_layers <- 3

results_file <- "results/predictions_island_scale_check.rds"

# read the prediction data if you already have it, and if not generate predictions

# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

if (file.exists(results_file)) {
  print("Existing results file found — reading the file and proceeding to analysis")
  
  combined_results <- readRDS(results_file)
  print("finished loading prediction results")
  
} else { # or alternatively run the prediction pipeline
  # Loop through all combinations of layers_to_train and layer_to_predict
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
      
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
        
        ### ---- b. + d. prediction with SVD and apply for all network combinations ----
        k_values <- c(2, 5, 10)
        lam0 <- lambda0(C)
        lambda_values <- c(1, 5, 50, 100, lam0)
        
        # Initialize variables to store the best results
        results <- data.frame(k = integer(),
                              lambda = numeric(),
                              original_links = numeric(),
                              predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL
        
        # Loop over all combinations of k and lambda
        for (k in k_values) {
          for (lambda in lambda_values) {
            # imputation
            r <- implement_impute(C, k, lambda)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
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
  # save the results
  saveRDS(combined_results, file = results_file)
} 

# after reading or producing the results, filter these (important!):
combined_results <- combined_results %>% 
  filter(k == 2) %>% 
  filter(!(input_lambda %in%  c(1, 5, 50, 100)))

## ---- 2. analysis ----
summary(combined_results)

combined_results_code <- combined_results

# claude version
implement_impute_P <- function(C, k, lambda) {
  # Apply softImpute
  
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  
  # Debias the fit to remove regularization effects
  # fit <- deBias(C, fit)
  
  # Reconstruct the matrix
  C_reconstructed <- softImpute::complete(C, fit)
  
  # Extract the reconstructed P matrix from C_reconstructed
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
  # Back-transform to original scale
  P_reconstructed_orig <- P_reconstructed +
    outer(row_centers[rownames(P)], col_centers[colnames(P)], "+")
  
  
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
  predicted_values <- P_reconstructed_orig[cbind(test_rows, test_cols)]
  
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

set.seed(42)
results_file <- "results/predictions_island_scale_check_P.rds"

# read the prediction data if you already have it, and if not generate predictions

# Initialize a data frame to store combined results for all layer combinations
combined_results <- data.frame()

if (file.exists(results_file)) {
  print("Existing results file found — reading the file and proceeding to analysis")
  
  combined_results <- readRDS(results_file)
  print("finished loading prediction results")
  
} else { # or alternatively run the prediction pipeline
  # Loop through all combinations of layers_to_train and layer_to_predict
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
      
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
        
        ### ---- b. + d. prediction with SVD and apply for all network combinations ----
        k_values <- c(2, 5, 10)
        lam0 <- lambda0(C)
        lambda_values <- c(1, 5, 50, 100, lam0)
        
        # Initialize variables to store the best results
        results <- data.frame(k = integer(),
                              lambda = numeric(),
                              original_links = numeric(),
                              predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL
        
        # Loop over all combinations of k and lambda
        for (k in k_values) {
          for (lambda in lambda_values) {
            # imputation
            r <- implement_impute_P(C, k, lambda)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
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
  # save the results
  saveRDS(combined_results, file = results_file)
} 

# after reading or producing the results, filter these (important!):
combined_results <- combined_results %>% 
  filter(k == 2) %>% 
  filter(!(input_lambda %in%  c(1, 5, 50, 100)))

## ---- 2. analysis ----
summary(combined_results)

combined_results_P <- combined_results
identical(combined_results_P, combined_results_code)
