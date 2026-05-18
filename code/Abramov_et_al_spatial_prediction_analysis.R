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
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
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
      axis.text.x     = element_text(size = 14)                          )
  
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
num_layers <- length(unique(aggregated_df$layer_from))

results_file <- "results/predictions_island_scale.rds"

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

# convert negatives to zeros
df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# know thy network - what species do we have in the system?
plant_species <- unique(df$node_from)           # get unique names
pollinator_species <- unique(df$node_to)
length(plant_species)
length(pollinator_species)

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

df_avg_plot <- df_avg %>% mutate(metric = recode(metric,
                                                 specificity       = "Specificity",
                                                 precision         = "Precision",
                                                 recall            = "Recall",
                                                 f05_score         = "F0.5 score",
                                                 balanced_accuracy = "Balanced accuracy",
                                                 mcc               = "MCC"
))

# 4) plot
optimal_threshold <- ggplot(df_avg_plot, aes(threshold, value, color = metric)) +
  geom_line(size = 1) +
  labs(
    x     = "Probability threshold",
    y     = "Average metric",
    color = "Metric"
  ) +
  scale_color_brewer(palette = "Pastel2") +
  tme

# pdf(
#   file   = "results/paper_figs/optimal_threshold.pdf",
#   width  = 5,    # inches
#   height = 4,
#   family = "Helvetica"   # or another installed font
# )
# print(optimal_threshold)
# dev.off()     # close the file

df_eval <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1L, 0L)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    # ROC-AUC (coerce to numeric!)
    auc_roc = tryCatch({
      roc_obj <- roc(response = original_binary,
                     predictor = predicted_prob,
                     quiet = TRUE, na.rm = TRUE,
                     levels = c(0,1), direction = "<")
      as.numeric(auc(roc_obj))   # <-- important
    }, error = function(e) NA_real_),
    
    # PR-AUC (guard against all-one-class cases)
    auc_pr = tryCatch({
      pos <- predicted_prob[original_binary == 1]
      neg <- predicted_prob[original_binary == 0]
      if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
      pr_obj <- pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)
      pr_obj$auc.integral
    }, error = function(e) NA_real_)
  ) %>%
  ungroup()

df_eval_summary <- df_eval %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    auc_roc_mean = mean(auc_roc, na.rm = TRUE),
    auc_roc_sd   = sd(auc_roc,   na.rm = TRUE),
    auc_pr_mean  = mean(auc_pr,  na.rm = TRUE),
    auc_pr_sd    = sd(auc_pr,    na.rm = TRUE),
    .groups = "drop"
  )


# 1) pivot to wide so F0.5 and balanced_accuracy are columns
df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# 2) find the threshold with the optimal f05 score
best_discrete <- df_wide %>%
  slice_max(f05_score, n = 1)

# results:
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold

df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

### ---- Fig. S13: predicted vs. observed weights ----
predicted_original <- df_removed %>%
  ggplot(aes(x = original_links, y = predicted_values)) +
  geom_point(alpha = 0.6, color = "lightsteelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linetype = "dashed") +
  stat_cor(method = "pearson", label.x = 55, label.y = 50) +  # change method to "spearman" if needed
  labs(
    x = "Weight of original links",
    y = "Predicted values"
    #title = "Correlation between predictions and original links"
  ) +
  geom_abline (slope=1, linetype = "dashed", color="salmon")+
  coord_equal()+
  theme_minimal(base_size = 11) + tme

pdf(
  file   = "results/paper_figs/predicted_original.pdf",
  width  = 6,    # inches
  height = 9,
  family = "Helvetica"   # or another installed font
)
print(predicted_original)
dev.off()     # close the file

png(
  filename = "results/paper_figs/predicted_original.png",
  width    = 6,     # inches
  height   = 9,
  units    = "in",
  res      = 300,   # resolution (dpi)
  type     = "cairo"  # better text rendering (recommended)
)

print(predicted_original)

dev.off()

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

### ---- Fig. S11: non-thresholded evaluation ----
# first add layer names
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# create a new grouped tibble for island names
new_layer_names <- net_name %>%
  mutate(group_id = (layer_id + 1) %/% 2) %>%  # Group pairs into 1, 2, 3...
  group_by(group_id) %>%
  summarise(name = gsub(" site.*", "", first(name)), .groups = "drop")  # Keep only location name

# Add names to main table
df_eval_summary <- df_eval_summary %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)

# now heatmaps

# roc
lims_roc <- range(df_eval_summary$auc_roc_mean, na.rm = TRUE)
mid_val_roc <- mean(lims_roc)

island_heatmap_auc <- 
  ggplot(df_eval_summary, aes(x = train_layer_name, y = test_layer_name, fill = auc_roc_mean)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = df_eval_summary[df_eval_summary$train_layer == df_eval_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "rosybrown2", 
                       midpoint = mid_val_roc, na.value = "gray") +  # Set NA values to gray
  labs(x = "Auxiliary location", y = "Target location", fill = "ROC-AUC") +
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

print(island_heatmap_auc)

# pr
lims_pr <- range(df_eval_summary$auc_pr_mean, na.rm = TRUE)
mid_val_pr <- mean(lims_pr)

island_heatmap_pr <- 
  ggplot(df_eval_summary, aes(x = train_layer_name, y = test_layer_name, fill = auc_pr_mean)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = df_eval_summary[df_eval_summary$train_layer == df_eval_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "thistle", 
                       midpoint = mid_val_pr, na.value = "gray") +  # Set NA values to gray
  labs(x = "Auxiliary location", y = "Target location", fill = "PR-AUC") +
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

print(island_heatmap_pr)

# combine 

pr_roc <- plot_grid(
  island_heatmap_auc + theme(plot.margin = unit(c(0,0,0,0), "cm")),
  island_heatmap_pr + theme(axis.title.y = element_blank()),
  rel_widths = c(1, 0.92),
  labels = c("(a)", "(b)"),
  label_size = 18,
  label_y = 0.8  # adjust this (e.g., 0.95, 0.9) to move labels closer
)

# supplementary figure pr_roc
pdf(file   = "results/paper_figs/pr_roc.pdf",
    width  = 13,    # inches
    height = 10,
    family = "Helvetica"   # or another installed font
)
pr_roc
dev.off()

### ---- Fig. S12: false positive rate ----
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
    title = "ROC curve",
    subtitle = paste0("AUC = ", round(auc_value, 3)),
    x = "False positive rate (1 − specificity)",
    y = "True positive rate (sensitivity)"
  ) +
  coord_equal() +
  theme_classic(base_size = 14) + tme

pdf(file   = "results/paper_figs/roc_curve.pdf",
    width  = 7,    # inches
    height = 7,
    family = "Helvetica"   # or another installed font
)
roc_curve
dev.off()

### ---- Fig. 2d: distribution of evaluators with/without external data ----
# this analysis shows us if predictions made using added information from other locations (off-diagonals in layer-to-layer predictions, as a heatmap) is any better than not adding any information (cases on the diagonal)

result_summary <- result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Single location",
    train_layer != test_layer ~ "Added location"
  ))

custom_colors <- c("Single location" = "steelblue",
                   "Added location" = "thistle")

# plot the histogram: 


hist_f05a <- plot_hist(
  result_summary,
  metric = "f05_score",
  y_axis_label = "Count of instances",
  x_axis_label = "F0.5 score"
) +
  scale_fill_manual(
    values = custom_colors,
    labels = c(
      "Single location" = "Within location",
      "Added location" = "External location"
    )
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.02)) +
  theme(
    axis.text.x = element_text(size = 22, hjust = 0.5),
    axis.text.y = element_text(size = 22),
    axis.title  = element_text(size = 30),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 21)
  )

hist_f05a # Fig. 2b

results_diags <- result_summary %>% filter(layer_comparison == "Single location")
results_offs <- result_summary %>% filter(layer_comparison == "Added location")

range(results_diags$f05_score)
mean(results_diags$f05_score)
sd(results_diags$f05_score)

range(results_offs$f05_score)
mean(results_offs$f05_score)
sd(results_offs$f05_score)

# # Base‐R PDF device
pdf(
  file   = "results/paper_figs/hist_f05a_legend_bottom.pdf",
  width  = 7.5,    # inches
  height = 7,
  family = "Helvetica"   # or another installed font
)
print(hist_f05a)
dev.off()     # close the file

# stats

# do we need welch/wilcoxon?
# first normality check
result_summary %>%
  group_by(layer_comparison) %>%
  shapiro_test(f05_score) # distributions are normal

# variance check
result_summary %>% levene_test(f05_score ~ layer_comparison)
# variances are similar. but due to the dependency between observations (not paired) we will not use Wilcoxon's or t-test

### ---- e. ecological inference ----
### ---- network size and density correlation with evaluators ----
# here we calculate the size and density of the networks and correlate them with evaluation metrics.
#### ---- calculate the size and density of our networks ----
# Initialize a data frame to store combined results for all layer combinations
results <- data.frame()

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {
    
    print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
    
    # Build the aggregated matrix A for training
    A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
    
    # Build the layer to predict matrix P
    P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
    
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

# summerize (Table S2)

# Add island names to main table
result_summary <- result_summary %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)

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
overall_mean_density <- mean(df_summary$density_P, na.rm = TRUE)

#### ---- Fig. S19: correlate network size with evaluators ----

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
  file   = "results/paper_figs/netsize_f05_nnse.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(netsize_f05_nnse)
dev.off()     # close the file

#### ---- Fig. S18: density ----

df_f05_nnse_density <- result_summary %>%
  select(f05_score, nnse, density_P, density_C) %>%
  pivot_longer(cols = c(density_P, density_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f05_score, nnse), names_to = "evaluator", values_to = "evaluator_value")

netdensity_f05_nnse <- plot_f05_nnse_vs_density_free_both(df_f05_nnse_density) + tme
netdensity_f05_nnse + theme(axis.text.x = element_text(size = 12))
netdensity_f05_nnse

pdf(
  file   = "results/paper_figs/netdensity_f05_nnse.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(netdensity_f05_nnse)
dev.off()     # close the file

### ---- Jaccard correlation with evaluators ----
#### ---- calculate Jaccard ----
# here we calculate the Jaccard index for every pair of networks
# Initialize a data frame to store combined results for all layer combinations
results_jaccard <- data.frame()

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

#### ---- Fig. 4a,b,c: plot Jaccard ----
# we filter only pairs of different islands for this analysis, since Jaccard index for the same island is 1 for all islands
canary_results_jaccard <- result_summary %>%
  filter(train_layer != test_layer)

jaccard_isl_f05 <- make_facet_scatter_plot(data = canary_results_jaccard, 
                                          evaluator = "f05_score",
                                          pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                          x_lab = "Jaccard similarity",
                                          y_lab = "F0.5 score",
                                          facet_scales = "free_x")
jaccard_isl_f05 # Fig. 4

# # Base‐R PDF device
# pdf(
#   file   = "jaccard_isl_f05.pdf",
#   width  = 7,    # inches
#   height = 3.5,
#   family = "Helvetica"   # or another installed font
# )
# print(jaccard_isl_f05)
# dev.off()     # close the file

jaccard_isl_nnse <- make_facet_scatter_plot(data = canary_results_jaccard, 
                                            evaluator = "nnse",
                                            pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                            x_lab = "Jaccard similarity",
                                            y_lab = "NNSE",
                                            facet_scales = "free_x")
jaccard_isl_nnse

### ---- degree impact on link assignment ----
# here we examine if the algorithm assigns more links to species with higher degree.

#### ---- calculate overall degree per species ----
# filter only existing links
df_filtered <- df %>%
  filter(itr == 1, original_links != 0) # filter existing interactions

# calculate degree for each plant species (node_from)
plant_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(plant_degree = n(), .groups = "drop")

# average plant degree by train_layer and test_layer
avg_plant_degree <- plant_degree %>%
  group_by(node_from) %>%
  summarise(avg_plant_degree = mean(plant_degree), .groups = "drop")

overall_plant_degree <- df_filtered %>% 
  group_by(node_from) %>% 
  summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")

# calculate degree for each pollinator species (node_to)
pollinator_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(poll_degree = n(), .groups = "drop")

# average pollinator degree by train_layer and test_layer
avg_pollinator_degree <- pollinator_degree %>%
  group_by(node_to) %>%
  summarise(avg_pollinator_degree = mean(poll_degree), .groups = "drop")

overall_poll_degree <- df_filtered %>% 
  group_by(node_to) %>% 
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

#### ---- Fig. 3c: plot degree vs. number of never observed interactions ----
# here by "island" we refer to a layer pair

df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

# for each island and interaction, determine if the interaction was observed.
# we use `any(original_links == 1)` so that if the interaction is observed in at least one iteration, we count it.
df_island <- df %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )

table(df_island$observed)

# now, for each unique interaction, compute:
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

# filter the interactions that were never observed throughout the data set
df_never_observed <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop") #plants

df_never_observed_poll <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")

# for relating to overall degree

never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")

never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")

# calculate correlation
df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants

# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 2)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

# specify species you want to label, if any
from_label <- c("Euphorbia_balsamifera_m",
                "Euphorbia_balsamifera_f",
                "Launaea_arborescens")

# build a little data‐frame just for the labels
df_labels <- df_to_correlate %>%
  filter(node_from %in% from_label) %>%
  distinct(node_from, .keep_all = TRUE) %>%
  mutate(
    label = paste0(
      'italic("',
      gsub("_", " ", node_from),
      '")'
    )
  )

# plot

plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "#2E8B57") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, \nnon-observed interactions",
    title = paste("Plants:", label_text_plants)   # <--- add label in title
  ) +
  theme_minimal() + tme +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# add repel‐text layer
plant_degree <- plant_degree +
  geom_text_repel(
    data      = df_labels,
    aes(label = label),
    parse     = TRUE,       # interpret the label as an expression
    size      = 4,          # tweak text size as needed
    box.padding   = 0.35,   # how much to push labels away from each other
    point.padding = 0.5,    # how much to push labels away from the points
    nudge_y       = -0.2     # optional small shift upward
  )

# same for pollinators
df_to_correlate_poll <- never_poll_degree_overall
df_to_correlate_poll$x <- df_to_correlate_poll$overall_poll_degree
df_to_correlate_poll$y <- df_to_correlate_poll$count_never_observed

correlation_poll <- cor.test(df_to_correlate_poll$x, df_to_correlate_poll$y, use = "complete.obs", method = "pearson")
correlation_poll

# Extract correlation coefficient and p-value
r_value <- round(correlation_poll$estimate, 2)
p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_to_correlate_poll, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "#CC8844") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, \nnon-observed interactions",
    title = paste("Pollinators:", label_text_polls)   # <--- add label in title
  ) +
  theme_minimal() + tme +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# # Create the figure
# final_plot <- combine_plots(plant_degree, poll_degree) # Fig. 3c

# Your two plots (remove individual axis labels)
plant_degree_clean <- plant_degree +
  labs(x = NULL, y = NULL)

poll_degree_clean <- poll_degree +
  labs(x = NULL, y = NULL)

# Add bottom label with padding
final_plot <- plot_grid(
  # main plots
  plot_grid(plant_degree_clean, poll_degree_clean, ncol = 2, align = "hv"),
  # x label
  ggdraw() + draw_label("Overall degree", fontface = "bold", size = 16),
  ncol = 1,
  rel_heights = c(1, 0.08)  # second element is space for x-axis label
)

# Add y label with padding
final_plot <- plot_grid(
  ggdraw() + draw_label("Number of predicted,\nnon-observed links",
                        angle = 90, fontface = "bold", size = 16),
  final_plot,
  ncol = 2,
  rel_widths = c(0.08, 1)   # first element is space for y-axis label
)

final_plot # fig. 3c

# # Save to PDF
pdf("results/paper_figs/degree_unobserved_links.pdf", width = 10, height = 7)  # adjust size as needed
grid::grid.draw(final_plot)
dev.off()

### ---- Fig. 3a: mapping never-observed links ----
# here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of cases in which these interactions were observed.

# order species by their degree
df_summary <- df_summary %>% left_join(overall_poll_degree, by="node_to")
df_summary <- df_summary %>% left_join(overall_plant_degree, by="node_from")

# determine the order of species in the plot based on their degree
plant_order <- df_summary %>%
  distinct(node_from, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(node_from)

poll_order <- df_summary %>%
  distinct(node_to, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(node_to)

# calculate proportion of islands in which each interaction occurs, rather than island pairs
df_filtered_self <- df_filtered %>% filter(train_layer == test_layer)

result_count <- df_filtered_self %>%
  # group by the interaction
  group_by(node_to, node_from) %>%
  summarise(
    # count unique test layers for this interaction
    n_test_layers = n_distinct(test_layer),
    .groups = "drop"
  ) %>%
  # calculate the proportion
  mutate(
    prop_test_layers = n_test_layers / n_distinct(df$test_layer)
  )

print(result_count)

# join to df_summary

df_summary <- df_summary %>%
  left_join(result_count %>% select(node_to, node_from, prop_test_layers),
            by = c("node_to", "node_from")) %>%
  # replace avg_prop with the calculated proportion
  mutate(avg_prop_isl = if_else(is.na(prop_test_layers), 0, prop_test_layers)) %>%
  select(-prop_test_layers)  # remove helper column if not needed

all((df_summary$avg_prop == 0) == (df_summary$avg_prop_isl == 0)) # check

# reset the levels for the species factors
df_summary$node_from <- factor(df_summary$node_from, levels = plant_order)
df_summary$node_to   <- factor(df_summary$node_to, levels = poll_order)

# map_missing_links <- ggplot(df_summary, aes(x = node_to, y = node_from)) +
#   # First layer: background heatmap for proportion observed (blue gradient)
#   geom_tile(aes(fill = avg_prop_isl)) +
#   scale_fill_gradient(low = "white", high = "steelblue", 
#                       name = "Observed links:\nproportion\nof islands\nobserved",
#                       breaks = seq(0, 1, 0.2)) +
#   
#   # Reset fill scale so the next layer can have its own gradient
#   new_scale_fill() +
#   
#   # Second layer: overlay only cells that were never observed but have high predicted value
#   geom_tile(
#     data = df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
#     aes(fill = avg_sigm_predicted),
#     alpha = 0.6
#   ) +
#   scale_fill_gradient(low = "tan1", high = "tomato2", 
#                       name = "Predicted links:\naverage predicted\nprobability",
#                       breaks = seq(0, 1, 0.1)) +
#   
#   # Final adjustments
#   theme_minimal() +
#   labs(x = "Pollinator", y = "Plant") +
#   theme(
#     axis.text.x = element_blank(), 
#     axis.text.y = element_text(size = 10),
#     legend.text = element_text(size = 12),
#     legend.position = "bottom",         # Place legends at the bottom
#     legend.box = "horizontal" 
#   ) + tme +
#   scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
#     bquote(italic(.(paste(y, collapse = " "))))
#   }))
# 
# print(map_missing_links)

# pdf(
#   file   = "results/paper_figs/map_missing_links.pdf",
#   width  = 11,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(map_missing_links)
# dev.off()     # close the file

### ---- difference in links predicted with/without external data ----
# this analysis shows us which links (and how many) were predicted only using external data, single-island data or combination of both.
df_island_sep <- df_island %>%
  separate(island_id, into = c("island1", "island2"), sep = "_", convert = TRUE)

df_island_diag <- df_island_sep %>%
  filter(island1 == island2)

# only within-island data
df_summary_diag <- df_island_diag %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

# only with extrnal data
df_island_offs <- df_island_sep %>%
  filter(island1 != island2)

df_summary_offs <- df_island_offs %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

identical(df_summary_offs$node_from, df_summary_diag$node_from) # check

# create a combined table
diff_df <- full_join(
  df_summary_offs,
  df_summary_diag,
  by     = c("node_from","node_to"),
  suffix = c("_offs","_diag")
) 

# these are interactions that were predicted to exist in the off-diagonals but not in the diagonal
diff_df_filtered <- diff_df %>% filter(diff_df$avg_sigm_predicted_diag < best_discrete_threshold & diff_df$avg_sigm_predicted_offs >= best_discrete_threshold)

# and these are interactions that were predicted to exist in the diagonal but not in the off-diagonals
diff_df_diags <- diff_df %>% filter(avg_sigm_predicted_offs < best_discrete_threshold & avg_sigm_predicted_diag >= best_discrete_threshold)

# order species by their degree
diff_df$node_from <- factor(diff_df$node_from, levels = plant_order)
diff_df$node_to   <- factor(diff_df$node_to, levels = poll_order)

# plot the differences
df_plot <- diff_df %>%
  filter(avg_prop_diag == 0) %>% # non-observed interactions
  mutate(
    sigm_cat = case_when(
      avg_sigm_predicted_diag < best_discrete_threshold &
        avg_sigm_predicted_offs >= best_discrete_threshold ~ "offs↑ only",
      
      avg_sigm_predicted_offs < best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "diag↑ only",
      
      avg_sigm_predicted_offs >= best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "both↑",
      
      TRUE ~ NA_character_
    )
  )

# map_missing_links_diags_offs <- ggplot(df_plot, aes(x = node_to, y = node_from)) +
#   
#   # allow a second fill scale
#   new_scale_fill() +
#   geom_tile(
#     data  = filter(df_plot, !is.na(sigm_cat)),
#     aes(fill = sigm_cat),
#     alpha = 0.6
#   ) +
#   scale_fill_manual(
#     values = c(
#       "offs↑ only" = "salmon",
#       "diag↑ only" = "plum3",
#       "both↑"       = "lightsteelblue"
#     ),
#     na.value = NA,
#     name   = "Difference in \nprediction approach",
#     labels = c(
#       "offs↑ only" = "Predicted only by \nadding external location",
#       "diag↑ only" = "Predicted only by \nsingle location",
#       "both↑"       = "Predicted by both approaches"
#     )
#   ) +
#   
#   # tidy up
#   theme_minimal() +
#   labs(x = "Pollinator", y = "Plant") +
#   theme(
#     axis.text.x  = element_blank(),
#     axis.text.y  = element_text(size = 8),
#     legend.position = "bottom",
#     legend.box      = "vertical"
#   ) +
#   
#   # your italic‐species labels and extra theme element
#   scale_y_discrete(
#     labels = function(x) lapply(strsplit(x, "_"), function(y) {
#       bquote(italic(.(paste(y, collapse = " "))))
#     })
#   ) +
#   tme
# 
# map_missing_links_diags_offs
# 
# pdf(
#   file   = "results/paper_figs/map_missing_links_diags_offs.pdf",
#   width  = 11,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(map_missing_links_diags_offs)
# dev.off()     # close the file

# create a data frame for the plot (unifies previous Fig. 3a and Fig. S11)
df_combined <- df_summary %>%
  select(node_from, node_to, avg_prop_isl) %>%
  full_join(
    df_plot %>%
      select(
        node_from, node_to,
        avg_sigm_predicted_offs,
        avg_sigm_predicted_diag,
        sigm_cat
      ),
    by = c("node_from", "node_to")
  ) %>%
  mutate(
    observed_link = avg_prop_isl > 0,
    
    pred_prob = case_when(
      sigm_cat == "offs↑ only" ~ avg_sigm_predicted_offs,
      sigm_cat == "diag↑ only" ~ avg_sigm_predicted_diag,
      sigm_cat == "both↑" ~ rowMeans(cbind(avg_sigm_predicted_offs, avg_sigm_predicted_diag), na.rm = TRUE), # if a link was predicted by both approaches, show the average probability
      TRUE ~ NA_real_
    ),
    
    sigm_cat_plot = case_when(
      avg_prop_isl == 0 & sigm_cat == "diag↑ only" ~ "Single location",
      avg_prop_isl == 0 & sigm_cat == "offs↑ only" ~ "External location",
      avg_prop_isl == 0 & sigm_cat == "both↑" ~ "Both approaches",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    sigm_cat_plot = factor(
      sigm_cat_plot,
      levels = c("Single location", "External location", "Both approaches")
    ),
    pred_prob_alpha = ifelse(
      !is.na(pred_prob),
      rescale(pmax(pred_prob, best_discrete_threshold),
              to = c(0.35, 1),
              from = c(best_discrete_threshold, 1)),
      NA_real_
    )
  )

# saveRDS(df_combined, "results/df_combined.rds")

n_unobserved_predicted <- df_combined %>%
  filter(
    avg_prop_isl == 0,
    !is.na(pred_prob),
    pred_prob > best_discrete_threshold
  ) %>%
  nrow()


# plot 

map_missing_links_merged <- ggplot(df_combined, aes(x = node_to, y = node_from)) +
  
  # observed links
  geom_tile(
    data = df_combined %>% filter(observed_link),
    aes(fill = avg_prop_isl)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "mediumaquamarine",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    name = "Observed links:\nproportion of islands observed",
    guide = guide_colorbar(
      order = 1,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),
      barheight = unit(0.45, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.6
    )
  ) +
  
  # black borders around observed links
  geom_tile(
    data = df_combined %>% filter(observed_link),
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  
  # reset fill scale
  new_scale_fill() +
  
  # predicted missing links:
  # fill = category hue, alpha = predicted probability
  geom_tile(
    data = df_combined %>%
      filter(!is.na(sigm_cat_plot), !is.na(pred_prob), pred_prob >= best_discrete_threshold),
    aes(fill = sigm_cat_plot, alpha = pred_prob),
    color = NA
  ) +
  
  scale_fill_manual(
    values = c(
      "Single location" = "plum3",
      "External location" = "salmon",
      "Both approaches" = "lightsteelblue"
    ),
    labels = c(
      "Single location" = "Within location",
      "External location" = "External location",
      "Both approaches" = "Both approaches"
    ),
    name = "Prediction approach",
    guide = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(alpha = 1)
    )
  ) +
  
  scale_alpha_continuous(
    limits = c(best_discrete_threshold, 1),
    range = c(0.35, 1),
    breaks = seq(best_discrete_threshold, 1, 0.1),
    oob = squish,
    name = "Predicted probability",
    guide = guide_legend(
      order = 3,
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1
    )
  ) +
  
  labs(x = "Pollinator", y = "Plant") +
  theme_minimal() + tme +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 24)),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(margin = margin(t = 24)),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.spacing.x = unit(0.2, "cm"),
    legend.box.spacing = unit(0.05, "cm")
  ) +
  scale_y_discrete(
    labels = function(x) lapply(strsplit(x, "_"), function(y) {
      bquote(italic(.(paste(y, collapse = " "))))
    })
  ) 

map_missing_links_merged

pdf(
  file   = "results/paper_figs/map_missing_links_merged.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(map_missing_links_merged)
dev.off()     # close the file


# how many links did each category add?
df_plot %>%
  filter(avg_prop_diag == 0) %>% 
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  )

### ---- existing predicted interactions ----
# if we want to know how each category contributed to verified existing links (that were observed in the system)
# plot the differences
df_plot_verified <- diff_df %>%
  filter(avg_prop_diag != 0) %>% # observed interactions
  mutate(
    sigm_cat = case_when(
      avg_sigm_predicted_diag < best_discrete_threshold &
        avg_sigm_predicted_offs >= best_discrete_threshold ~ "offs↑ only",
      
      avg_sigm_predicted_offs < best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "diag↑ only",
      
      avg_sigm_predicted_offs >= best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "both↑",
      
      TRUE ~ NA_character_
    )
  )

# # if we want to plot them
# map_existing_links_predicted <- ggplot(df_plot_verified, aes(x = node_to, y = node_from)) +
# 
#   # allow a second fill scale
#   new_scale_fill() +
#   geom_tile(
#     data  = filter(df_plot_verified, !is.na(sigm_cat)),
#     aes(fill = sigm_cat),
#     alpha = 0.6
#   ) +
#   scale_fill_manual(
#     values = c(
#       "offs↑ only" = "salmon",
#       "diag↑ only" = "plum3",
#       "both↑"       = "lightsteelblue"
#     ),
#     na.value = NA,
#     name   = "Difference in \nprediction approach",
#     labels = c(
#       "offs↑ only" = "Predicted only by \nadding external location",
#       "diag↑ only" = "Predicted only by \nsingle location",
#       "both↑"       = "Predicted by both approaches"
#     )
#   ) +
# 
#   # tidy up
#   theme_minimal() +
#   labs(x = "Pollinator", y = "Plant") +
#   theme(
#     axis.text.x  = element_blank(),
#     axis.text.y  = element_text(size = 8),
#     legend.position = "bottom",
#     legend.box      = "vertical"
#   ) +
# 
#   # your italic‐species labels and extra theme element
#   scale_y_discrete(
#     labels = function(x) lapply(strsplit(x, "_"), function(y) {
#       bquote(italic(.(paste(y, collapse = " "))))
#     })
#   ) +
#   tme
# 
# map_existing_links_predicted

# how many links did each category add?
df_plot_verified %>%
  filter(avg_prop_diag != 0) %>%
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  ) # 178 observed links were not predicted

df_plot_verified %>%
  filter(avg_prop_diag != 0) %>%
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  ) %>%
  filter(!is.na(sigm_cat)) %>%  # Exclude NA category
  mutate(
    total_links = sum(n_links),  # Sum of non-NA categories only
    percentage = n_links / total_links * 100   # Calculate percentage
  )

#### ---- how many of the predicted links have evidence ----
# so overall we had
# predicted links that were observed at least once in the system
# (df_plot_verified)
df_observed <- df_plot_verified %>%
  filter(avg_prop_diag != 0) %>%
  group_by(sigm_cat) %>%
  summarise(
    n_observed = n()
  ) %>%
  mutate(status = "observed")  # Add a column to mark these as observed

# predicted links that were never observed in the system
df_non_observed <- df_plot %>%
  filter(avg_prop_diag == 0) %>%
  group_by(sigm_cat) %>%
  summarise(
    n_non_observed = n()
  ) %>%
  mutate(status = "non_observed")  # Add a column to mark these as non-observed

# Combine both dataframes
df_combined_obs_non <- bind_rows(df_observed, df_non_observed)

# in total, the proportion of predicted links that were observed in the system:
df_combined_obs_non %>%
  filter(!is.na(sigm_cat)) %>%
  summarise(
    total_observed = sum(n_observed, na.rm = TRUE),
    total_non_observed = sum(n_non_observed, na.rm = TRUE)
  ) %>%
  mutate(
    total_links = total_observed + total_non_observed,  # Total links (observed + non-observed)
    proportion_observed = total_observed / total_links * 100  # Proportion of observed links
  )


# pie chart
# 1. Count how many interactions fall into each category
df_counts <- df_plot %>%
  filter(sigm_cat != "NA") %>% 
  count(sigm_cat, name = "n") %>%
  arrange(desc(sigm_cat)) %>%          # optional: control legend/order
  mutate(
    frac = n / sum(n),                 # fraction of total
    pct  = percent(frac)               # human‑readable percent
  )

# 2. Define your colors
my_cols <- c(
  "offs↑ only"  = "salmon",
  "diag↑ only"  = "plum3",
  "both↑"       = "lightsteelblue"
)

new_labels <- c(
  "offs↑ only" = "Predicted only by external location" ,
  "diag↑ only" = "Predicted only by single location" ,
  "both↑"      = "Predicted by both approaches"
)

### ---- Fig. 3 (complete): mapping missing links and updated pie chart ----
#### ---- Fig. 3b: pie chart ----

# add percentages
df_counts <- df_counts |>
  dplyr::mutate(
    sigm_cat = dplyr::recode(
      sigm_cat,
      "offs↑ only" = "External location",
      "diag↑ only" = "Within location",
      "both↑"      = "Both"
    )
  )

values <- df_counts$n

labels <- paste0(df_counts$sigm_cat,
                 "\n",
                 df_counts$n,
                 " (", df_counts$pct, ")")

pdf(
  file   = "results/paper_figs/pie_chart.pdf",
  width  = 6.5,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
# Draw pie without labels
pie(values,
    labels = NA,
    col = my_cols,
    border = "white")

# # Add labels
fractions <- values / sum(values)
cum_fractions <- cumsum(fractions)
mid_angles <- 1.6 * pi * (cum_fractions - fractions / 2)

label_radius <- rep(1.3, length(values))
i <- which(df_counts$sigm_cat == "Both")
label_radius[i] <- 1.4   # move this one further out
j <- which(df_counts$sigm_cat == "Within location")
label_radius[j] <- 0.9
k <- which(df_counts$sigm_cat == "External location")
label_radius[k] <- 0.82

text_colors <- c("salmon", "plum3", "lightsteelblue")

text(label_radius * cos(mid_angles),
     label_radius * sin(mid_angles),
     labels = labels,
     adj = ifelse(cos(mid_angles) > 0, 0, 1),
     col = text_colors,
     cex = 1.2,
     xpd = TRUE)

dev.off()     # close the file


# Read the pie chart PDF as an image (first page)
img <- magick::image_read_pdf("results/paper_figs/pie_chart.pdf", density = 300)

image_info(img)

img_cropped <- image_crop(img, geometry = "1950x1800+0+280") # "WIDTHxHEIGHT+LEFT+TOP"
new_width  <- 1950 - 0 # how much to crop from right
new_height <- 1800 - 550 # how much to crop from bottom

img_cropped <- image_crop(img_cropped, 
                          geometry = paste0(new_width, "x", new_height, "+0+0"))


# Convert to grob
pie_grob <- cowplot::ggdraw() + cowplot::draw_image(img_cropped)

bottom_row <- plot_grid(
  pie_grob,
  final_plot,
  rel_widths = c(0.9, 1.2),
  labels = c("(b)", "(c)"),
  label_size = 15
)


fig3 <- plot_grid(map_missing_links_merged + theme(plot.margin = unit(c(0.8,0.2,0.2,0.2), "cm")), 
                  bottom_row, labels = c('(a)', ''), 
                  ncol = 1, rel_heights = c(1.3, 0.7), label_size = 15)


pdf(file   = "results/paper_figs/missing_interactions_degree2.pdf",
    width  = 13,    # inches
    height = 12,
    family = "Helvetica"   # or another installed font
)
fig3
dev.off()

### ---- distance decay ----
#### ---- add distances and location names ----
distance_table <- read.csv("data/distance_between_sites_canary.csv", row.names = NULL)

# create new table with averaged distances at the island level
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

# print result
print(distance_island_table)

# modify the 'from' and 'to' columns in distance_island_table
distance_island_table <- distance_island_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

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

#### ---- Fig. 4d:  distance correlation with evaluators ----

# remove sites form within the same island - only use information from different islands for distance decay
result_summary_island_dif <- result_summary_island %>% filter(train_layer != test_layer)

# plot
cor_plot_dif_isl_f05  <- make_cor_plot(result_summary_island_dif, evaluator = "f05_score", extra_theme = tme) +
  theme(
    axis.text.x = element_text(size = 14)  # Adjust x-axis text size
  )

## a check for relationship with log distance
make_cor_plot_log <- function(data, evaluator, 
                              distance_col = "log_distance_km", 
                              x_lab = "Log geographic distance (km)",
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

result_summary_island_dif <- result_summary_island_dif %>% mutate(log_distance_km = log(distance_km))

cor_plot_dif_isl_f05_log  <- make_cor_plot_log(result_summary_island_dif, evaluator = "f05_score", extra_theme = tme) +
  theme(
    axis.text.x = element_text(size = 14)  # Adjust x-axis text size
  )
# the correlation is nearly identical.

# pdf(
#   file   = "cor_plot_dif_isl_f05.pdf",
#   width  = 4,
#   height = 4,
#   family = "Helvetica"
# )
# print(cor_plot_dif_isl_f05)
# dev.off()     # close the file


# NNSE and distance
cor_plot_dif_isl_nnse  <- make_cor_plot(result_summary_island_dif, evaluator = "nnse", extra_theme = tme)
saveRDS(result_summary_island_dif, "results/result_summary_island_dif.rds")

### ---- Table 1: MRM test and variable importance ----
result_summary_island_dif <- readRDS("results/result_summary_island_dif.rds")

# build square matrices of f0.5 and distance
#    (layers must be in the same order for rows & cols)

# a) get a list of all unique layers
layers <- sort(unique(c(result_summary_island_dif$train_layer, result_summary_island_dif$test_layer)))

# b) initialize empty matrices
f05_mat      <- matrix(NA, nrow=length(layers), ncol=length(layers),
                      dimnames=list(layers, layers))
dist_mat_km <- f05_mat

# c) fill in each cell [i,j] with the corresponding f05_score and distance_km
for(i in layers) for(j in layers) {
  # subset rows where train=i and test=j
  sub <- result_summary_island_dif[result_summary_island_dif$train_layer==i & result_summary_island_dif$test_layer==j, ]
  if(nrow(sub)==1) {
    f05_mat[i,j]      <- sub$f05_score
    dist_mat_km[i,j] <- sub$distance_km
  }
}

# d) because MRM uses symmetric distance matrices, average [i,j] & [j,i]
f05_sym      <- sym_average(f05_mat)
dist_sym_km <- sym_average(dist_mat_km)

# save matrices for later use
write.csv(f05_sym, "results/f05_distance_matrix.csv", row.names = TRUE)
write.csv(dist_sym_km, "results/distance_matrix_km.csv", row.names = TRUE)


# e) convert to “dist” objects (lower triangle)
dist_f05      <- as.dist(f05_sym)
dist_km      <- as.dist(dist_sym_km)

# 4. run the MRM
#    — this will regress the F0.5‐distance matrix on the geographic–distance matrix
set.seed(42)   # for reproducibility of permutations
mrm_out <- MRM(dist_f05 ~ dist_km, nperm=999)

# 5. results
print(mrm_out)


# 6. now with composed MRM -> add C connectance, C matrix size, and Jaccard.
jaccard_mat <- matrix(NA, nrow=length(layers), ncol=length(layers),
                      dimnames=list(layers, layers))
c_connectance_mat <- jaccard_mat
c_mat_size_mat <- jaccard_mat

for(i in layers) for(j in layers) {
  sub <- result_summary_island_dif[result_summary_island_dif$train_layer==i & result_summary_island_dif$test_layer==j, ]
  if(nrow(sub)==1) {
    jaccard_mat[i,j]       <- sub$jaccard_edges
    c_connectance_mat[i,j] <- sub$density_C
    c_mat_size_mat[i,j]    <- sub$size_C
  }
}

# d) because MRM uses symmetric distance matrices, average [i,j] & [j,i]
jaccard_sym           <- sym_average(jaccard_mat)
c_connectance_sym     <- sym_average(c_connectance_mat)
c_mat_size_sym        <- sym_average(c_mat_size_mat)

# e) convert to “dist” objects (lower triangle)
jaccard           <- as.dist(jaccard_sym)
c_connectance     <- as.dist(c_connectance_sym)
c_mat_size        <- as.dist(c_mat_size_sym)

# now run MRM for each possible subset of predictors to see which ones are significant on their own, and which ones remain significant when controlling for others
# get every possible combination of predictors
predictors <- c("dist_km", "c_connectance", "c_mat_size", "jaccard")
predictor_combinations <- unlist(lapply(1:length(predictors), function(n) {
  combn(predictors, n, simplify = FALSE)
}), recursive = FALSE)

# for each combination, add dist_km predictor, including only dist_km as a baseline
#predictor_combinations <- 
#  lapply(predictor_combinations, function(preds) {c("dist_km", preds)})
#predictor_combinations <- c(list("dist_km"), predictor_combinations) # add dist_km alone as a baseline

# run MRM for each combination
mrm_results <- lapply(predictor_combinations, function(preds) {
  formula <- as.formula(paste("dist_f05 ~ ", paste(preds, collapse = " + ")))
  res <- MRM(formula, nperm=999)
  list(predictors = preds, result = res)
})
# extract results into a data frame
mrm_summary <- do.call(rbind, lapply(mrm_results, function(x) {
  data.frame(
    predictors = paste(x$predictors, collapse = " + "),
    r_squared = x$result$r.squared[1],
    p_value = x$result$r.squared[2]
  )
}))
print(mrm_summary)

##### ---- evaluate the best variable combination using AICc ----
# (aka Akaike Information Criterion corrected for small sample sizes)

# Note: distance-matrix entries are not independent, so AIC/AICc here 
# should be treated as a model-comparison heuristic (useful for ranking), 
# not a fully classical likelihood-based IC

# get the number of observations (number of unique pairs)
n <- length(jaccard) # number of unique pairs (lower triangle of the matrix)
# calculate AICc for each model
mrm_summary$aicc <- sapply(mrm_results, function(x) {
  k <- length(x$predictors) # number of predictors
  r2 <- x$result$r.squared[1] # R-squared of the model
  aic <- n * log(1 - r2) + 2 * k
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  return(aicc)
})
# rank models by AICc
mrm_summary <- mrm_summary[order(mrm_summary$aicc), ]
print(mrm_summary)

# so according to this, the best combination with lowest AICc is:
# dist_km + jaccard


# try to extract feature importance
# Compare R² difference when dropping predictors - Relative contribution to explained variance
# for each feature in "predictors" calculate the R² difference when dropping it from the full model
pred_full <- paste(predictors, collapse = " + ")
full_r2 <- mrm_summary$r_squared[mrm_summary$predictors == pred_full]

feature_importance <- sapply(predictors, function(pred) {
  # find the model that includes all predictors except the current one
  other_preds <- setdiff(predictors, pred)
  str <- paste(other_preds, collapse = " + ")
  r2_no_pred <- mrm_summary$r_squared[mrm_summary$predictors == str]
  importance <- full_r2 - r2_no_pred
  return(importance)
})

feature_importance_df <- data.frame(
  predictor = predictors,
  importance = feature_importance
) %>%
  arrange(desc(importance)) # most important predictor is Jaccard

# Add the feature importance column to mrm_summary
mrm_summary$importance <- sapply(mrm_summary$predictors, function(preds) {
  # Split predictors (e.g., "dist_km + jaccard") into individual predictors
  predictor_list <- strsplit(preds, " \\+ ")[[1]]
  
  # Sum the importance of all the predictors in the model
  sum_importance <- sum(feature_importance_df$importance[feature_importance_df$predictor %in% predictor_list])
  
  return(sum_importance)
})

print(mrm_summary)

### ---- Fig. 2c heatmap ----
# Compute limits and midpoint dynamically
lims <- range(result_summary_island$f05_score, na.rm = TRUE)
mid_val <- mean(lims)

island_heatmap_f05 <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = f05_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "salmon2", 
                       midpoint = mid_val, na.value = "gray") +  # Set NA values to gray
  labs(x = "Auxiliary location", y = "Target location", fill = "F0.5 score") +
  theme_minimal() + tme +
  theme(
    text = element_text(size = 14),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  # Rotate x-axis labels by 45 degrees
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  ) +
  coord_fixed()

print(island_heatmap_f05)

pdf(
  file   = "results/paper_figs/island_heatmap_f05.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(island_heatmap_f05)
dev.off()     # close the file

### ---- additional stats ----
overall_sd_f05 <- sd(result_summary$f05_score, na.rm = TRUE)
overall_mean_f05 <- mean(result_summary$f05_score, na.rm = TRUE)
overall_max_f05 <- max(result_summary$f05_score, na.rm = TRUE)


# sahara predictions
results_sahara <- result_summary_island %>% filter(test_layer_name == "Western Sahara" & train_layer_name != "Western Sahara")
mean(results_sahara$f05_score)
sd(results_sahara$f05_score)

# is including external data better?
# run t-test via formula interface
wilcox_f05 <- wilcox.test(f05_score ~ layer_comparison,
                         data = result_summary_island,
                         exact = FALSE)  # turn off exact test for larger samples
wilcox_f05

# do we need welch/wilcoxon?
# first normality check
shapiro_f05 <- result_summary_island %>%
  group_by(layer_comparison) %>%
  shapiro_test(f05_score) # for added location, the distribution is not normal

# variance check
levene_f05 <- result_summary_island %>% levene_test(f05_score ~ layer_comparison)
# variances are equal, use wilcoxon

## ---- Fig. S14: species occurrence and degree ----

# Combine plant and pollinator columns into one column of species occurrences

df_occurrence <- df_filtered %>%
  filter(test_layer == train_layer)

# For plants (node_to)
plant_occurrence <- df_occurrence %>%
  group_by(node_from) %>%
  summarise(num_islands = n_distinct(test_layer)) %>%
  arrange(desc(num_islands))

# For pollinators (node_from)
pollinator_occurrence <- df_occurrence %>%
  group_by(node_to) %>%
  summarise(num_islands = n_distinct(test_layer)) %>%
  arrange(desc(num_islands))

# View results
print(plant_occurrence)
print(pollinator_occurrence)

# join with degree data
# for plants
plant_data <- left_join(overall_plant_degree, plant_occurrence, by = "node_from")

# For pollinators
pollinator_data <- left_join(overall_poll_degree, pollinator_occurrence, by = "node_to")

# correlate
# calculate correlation
df_to_correlate_occurrence <- plant_data
df_to_correlate_occurrence$x <- df_to_correlate_occurrence$overall_plant_degree
df_to_correlate_occurrence$y <- df_to_correlate_occurrence$num_islands

correlation_plants_occ <- cor.test(df_to_correlate_occurrence$x, df_to_correlate_occurrence$y, use = "complete.obs", method = "pearson")
correlation_plants_occ

# Extract correlation coefficient and p-value
r_value_occ <- round(correlation_plants_occ$estimate, 2)
p_value_occ <- formatC(correlation_plants_occ$p.value, digits = 2)  # or round as you prefer
label_text_plants_occ <- paste0("r = ", r_value_occ, ", p = ", p_value_occ)

# plot occurrs

plant_degree_occurrence <- ggplot(df_to_correlate_occurrence, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of islands occurrs",
    title = paste("Plants:", label_text_plants_occ)   # <--- add label in title
  ) +
  theme_minimal() + tme

# same for pollinators
df_to_correlate_poll_occ <- pollinator_data
df_to_correlate_poll_occ$x <- df_to_correlate_poll_occ$overall_poll_degree
df_to_correlate_poll_occ$y <- df_to_correlate_poll_occ$num_islands

correlation_poll_occ <- cor.test(df_to_correlate_poll_occ$x, df_to_correlate_poll_occ$y, use = "complete.obs", method = "pearson")
correlation_poll_occ

# Extract correlation coefficient and p-value
r_value_poll <- round(correlation_poll_occ$estimate, 2)
p_value_poll <- formatC(correlation_poll_occ$p.value, digits = 2)  # or round as you prefer
label_text_polls_occ <- paste0("r = ", r_value_poll, ", p = ", p_value_poll)

poll_degree_occurrence <- ggplot(df_to_correlate_poll_occ, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "rosybrown2") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of islands occurrs",
    title = paste("Pollinators:", label_text_polls_occ)   # <--- add label in title
  ) +
  theme_minimal() + tme

# # Create the final figure
# Your two plots (remove individual axis labels)
plant_degree_occ_clean <- plant_degree_occurrence +
  labs(x = NULL, y = NULL)

poll_degree_occ_clean <- poll_degree_occurrence +
  labs(x = NULL, y = NULL)

# Add bottom label with padding
final_plot_occ <- plot_grid(
  plot_grid(
    plant_degree_occ_clean,
    NULL,
    poll_degree_occ_clean,
    ncol = 3,
    align = "hv",
    rel_widths = c(1, 0.04, 1),
    labels = c("(a)", "", "(b)"),
    label_size = 16,
    label_fontface = "bold",
    label_x = -0.02,
    label_y = 0.99,
    hjust = 0,
    vjust = 1
  ),
  ggdraw() + draw_label("Overall degree", fontface = "bold", size = 16),
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Add y label with padding
final_plot_occ <- plot_grid(
  ggdraw() + draw_label(
    "Number of islands present",
    angle = 90,
    fontface = "bold",
    size = 16
  ),
  final_plot_occ,
  ncol = 2,
  rel_widths = c(0.08, 1)
)

final_plot_occ

# # Save to PDF
pdf("results/paper_figs/degree_occurrence.pdf", width = 8, height = 5)  # adjust size as needed
grid::grid.draw(final_plot_occ)
dev.off()


# ---- local vs. global degrees and degree binning ----
## ---- Fig. S15: correlate predicted and observed local degree ----
# if we want to see if the amount of added links (locally in the layer combination) is correlated with the local degree of the species:
df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

# plants
# average predicted probability across iterations
pred_avg <- df %>%
  group_by(train_layer, test_layer, node_from, node_to, removed) %>%
  summarise(
    mean_pred = mean(predicted_prob_sigm, na.rm = TRUE),
    original_links = first(original_links),  # constant across itr
    .groups = "drop"
  )

# predicted degree per species per layer combo
pred_degree <- pred_avg %>%
  filter(removed == 1,
         mean_pred > best_discrete_threshold) %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(
    predicted_degree = n_distinct(node_to),
    .groups = "drop"
  )

# observed local degree
obs_degree <- pred_avg %>%
  filter(original_links > 0) %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(
    observed_degree = n_distinct(node_to),
    .groups = "drop"
  )

# combine
degree_comparison <- obs_degree %>%
  left_join(pred_degree,
            by = c("train_layer", "test_layer", "node_from")) %>%
  mutate(
    predicted_degree = replace_na(predicted_degree, 0)
  )

# correlate
correlation_plants_d <- cor.test(degree_comparison$observed_degree, degree_comparison$predicted_degree, use = "complete.obs", method = "pearson")
correlation_plants_d

# Extract correlation coefficient and p-value
r_value_d <- round(correlation_plants_d$estimate, 2)
p_value_d <- formatC(correlation_plants_d$p.value, digits = 2)  # or round as you prefer
label_text_plants_d <- paste0("r = ", r_value_d, ", p = ", p_value_d)

# plot 
plant_degree_obs_pred <- ggplot(degree_comparison, aes(x = observed_degree, y = predicted_degree)) +
  geom_point(alpha = 0.5, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Observed degree",
    y = "Predicted degree",
    title = paste("Plants:", label_text_plants_d)   # <--- add label in title
  ) +
  theme_minimal() + tme

plant_degree_obs_pred

# pollinators

# predicted degree per pollinator per layer combo (pollinators = node_to)
pred_degree_poll <- pred_avg %>%
  filter(removed == 1,
         mean_pred > best_discrete_threshold) %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(
    predicted_degree = n_distinct(node_from),
    .groups = "drop"
  )

# observed local degree per pollinator (FULL network; pollinators = node_to)
obs_degree_poll <- pred_avg %>%
  filter(original_links > 0) %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(
    observed_degree = n_distinct(node_from),
    .groups = "drop"
  )

# combine
degree_comparison_poll <- obs_degree_poll %>%
  left_join(pred_degree_poll,
            by = c("train_layer", "test_layer", "node_to")) %>%
  mutate(
    predicted_degree = replace_na(predicted_degree, 0)
  )

# correlate
correlation_poll_d <- cor.test(
  degree_comparison_poll$observed_degree,
  degree_comparison_poll$predicted_degree,
  use = "complete.obs",
  method = "pearson"
)
correlation_poll_d

# Extract correlation coefficient and p-value
r_value_d <- round(correlation_poll_d$estimate, 2)
p_value_d <- formatC(correlation_poll_d$p.value, digits = 2)
label_text_poll_d <- paste0("r = ", r_value_d, ", p = ", p_value_d)

# plot
poll_degree_obs_pred <- ggplot(degree_comparison_poll, aes(x = observed_degree, y = predicted_degree)) +
  geom_point(alpha = 0.5, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Observed degree",
    y = "Predicted degree",
    title = paste("Pollinators:", label_text_poll_d)
  ) +
  theme_minimal() + tme

poll_degree_obs_pred

# Remove individual axis labels
plant_degree_clean_d <- plant_degree_obs_pred +
  labs(x = NULL, y = NULL)

poll_degree_clean_d <- poll_degree_obs_pred +
  labs(x = NULL, y = NULL)

# Combine the two panels
main_panel <- plot_grid(
  plant_degree_clean_d,
  NULL,
  poll_degree_clean_d,
  ncol = 3,
  align = "hv",
  rel_widths = c(1, 0.04, 1),
  labels = c("(a)", "", "(b)"),
  label_size = 16,
  label_fontface = "bold",
  label_x = -0.02,
  label_y = 0.99,
  hjust = 0,
  vjust = 1
)


# Add bottom (shared x-axis) label
with_x_label <- plot_grid(
  main_panel,
  ggdraw() + draw_label("Observed degree",
                        fontface = "bold",
                        size = 16),
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Add left (shared y-axis) label
final_plot_d <- plot_grid(
  ggdraw() + draw_label("Predicted degree (withheld links)",
                        angle = 90,
                        fontface = "bold",
                        size = 16),
  with_x_label,
  ncol = 2,
  rel_widths = c(0.08, 1)
)

final_plot_d

# # Save to PDF
pdf("results/paper_figs/local_degree_predicted_links.pdf", width = 10, height = 7)  # adjust size as needed
grid::grid.draw(final_plot_d)
dev.off()

## ---- Fig. S17: plant degree across islands ----
df_itr1_observed <- df %>% filter(itr == 1 & original_links > 0)

plant_island_degree <- df_itr1_observed %>%
  group_by(test_layer, node_from) %>%
  summarise(island_degree = n_distinct(node_to), .groups = "drop")

plant_global_degree <- df_itr1_observed %>%
  group_by(node_from) %>%
  summarise(global_degree = n_distinct(node_to), .groups = "drop")

plant_degree_combined <- plant_island_degree %>%
  left_join(plant_global_degree, by = "node_from")

# 1) reorder species the same way as in the plot
plant_degree_plot_df <- plant_degree_combined %>%
  mutate(node_from = fct_reorder(node_from, island_degree, .fun = median))

# 2) count unique test layers per species
layer_counts <- plant_degree_plot_df %>%
  group_by(node_from) %>%
  summarise(
    n_test_layers = n_distinct(test_layer),
    y_pos = max(island_degree, na.rm = TRUE) + 2,   # place label above each box
    .groups = "drop"
  )

# 3) plot
plant_island_degree <- plant_degree_plot_df %>%
  ggplot(aes(x = node_from, y = island_degree, fill = node_from)) +
  geom_boxplot() +
  geom_text(
    data = layer_counts,
    aes(x = node_from, y = y_pos, label = n_test_layers),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "Plant", y = "Island level degree") +
  tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

plant_island_degree

pdf(
  file   = "results/paper_figs/plant_island_degree.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(plant_island_degree)
dev.off()     # close the file

## ---- Fig. S16: performance by degree binning ----
# add degrees to prediction data frame
df2 <- df %>%
  left_join(obs_degree,
            by = c("train_layer", "test_layer", "node_from"))

# create degree binning
df2 <- df2 %>%
  mutate(degree_bin = ntile(observed_degree, 3))

# and label them
df2$degree_bin <- factor(df2$degree_bin,
                         levels = 1:3,
                         labels = c("specialists", "intermediate", "generalists"))

# evaluate
df3 <- df2 %>%
  mutate(predicted_binary = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0))

# for each iteration:
perf_iter <- df3 %>%
  filter(removed == 1) %>% 
  group_by(itr, degree_bin) %>%
  summarise(
    TP = sum(original_links > 0 & predicted_binary == 1),
    FN = sum(original_links > 0 & predicted_binary == 0),
    TN = sum(original_links == 0 & predicted_binary == 0),
    FP = sum(original_links == 0 & predicted_binary == 1),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    F05 = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    .groups = "drop"
  )

# average across iterations:
perf_summary <- perf_iter %>%
  group_by(degree_bin) %>%
  summarise(
    mean_F05 = mean(F05, na.rm = TRUE),
    sd_F05 = sd(F05, na.rm = TRUE)
  )

# check for normality
perf_iter %>%
  group_by(degree_bin) %>%
  shapiro_test(F05) # normal, we can use anova

p_aov <- aov(F05 ~ degree_bin, data = perf_iter)
# Extract p-value
p_val <- summary(p_aov)[[1]]["degree_bin", "Pr(>F)"]

degree_binning <- ggplot(perf_iter, aes(x = degree_bin, y = F05, fill = degree_bin)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_brewer(palette = "Pastel2") + 
  labs(
    x = "Species degree class",
    y = expression(F[0.5]),
    title = paste0("ANOVA p = ", signif(p_val, 3))
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 12)
  ) + tme

pdf("results/paper_figs/degree_binning.pdf", width = 6, height = 6)  # adjust size as needed
grid::grid.draw(degree_binning)
dev.off()

# Combine key plots into figures -----------

### ---- Fig. 2: ----

# Fig. 2a is p_f05 from subset_analysis.R
# Fig. 2b is p_nnse from subset_analysis.R
# Fig. 2c is island_heatmap_f05
# Fig. 2d is hist_f05a

# Step 1: Read the PDF files (assuming your PDFs are stored at specific file paths)
p_f05_img <- magick::image_read_pdf("results/paper_figs/island_subset_f05.pdf", density = 300)  # Read the p_f05 PDF
p_nnse_img <- magick::image_read_pdf("results/paper_figs/island_subset_nnse.pdf", density = 300)  # Read the p_nnse PDF
hist_f05a_img <- magick::image_read_pdf("results/paper_figs/hist_f05a_legend_bottom.pdf", density = 300)  # Read the hist_f05a PDF
island_heatmap_f05_img <- magick::image_read_pdf("results/paper_figs/island_heatmap_f05.pdf", density = 300)  # Read the island_heatmap_f05 PDF

# crop
image_info(island_heatmap_f05_img)
island_heatmap_f05_img_cropped <- image_crop(island_heatmap_f05_img, geometry = "1800x1800+0+200") # "WIDTHxHEIGHT+LEFT+TOP"
new_width  <- 1800 - 0 # how much to crop from right
new_height <- 1800 - 400 # how much to crop from bottom

island_heatmap_f05_img_cropped <- image_crop(island_heatmap_f05_img_cropped, 
                                             geometry = paste0(new_width, "x", new_height, "+0+0"))


# Step 2: Convert the PDFs to graphical objects (grobs)
p_f05_grob <- cowplot::ggdraw() + cowplot::draw_image(p_f05_img)
p_nnse_grob <- cowplot::ggdraw() + cowplot::draw_image(p_nnse_img)
hist_f05a_grob <- cowplot::ggdraw() + cowplot::draw_image(hist_f05a_img)
island_heatmap_f05_grob <- cowplot::ggdraw() + cowplot::draw_image(island_heatmap_f05_img_cropped)


# Step 3: Add labels manually using ggdraw and draw_label
p_f05_grob_labeled <- p_f05_grob + cowplot::draw_label("(a)", x = 0, y = 1, hjust = 0, vjust = 1, size = 18, fontface = "bold")
p_nnse_grob_labeled <- p_nnse_grob + cowplot::draw_label("(b)", x = 0, y = 1, hjust = 0, vjust = 1, size = 18, fontface = "bold")
island_heatmap_f05_grob_labeled <- island_heatmap_f05_grob + cowplot::draw_label("(c)", x = 0, y = 1.1, hjust = 0, vjust = 1, size = 18, fontface = "bold")
hist_f05a_grob_labeled <- hist_f05a_grob + cowplot::draw_label("(d)", x = 0, y = 1.1, hjust = 0, vjust = 1, size = 18, fontface = "bold")


# Step 4: Create the 2x2 grid with labeled plots
fig2x2_labeled <- plot_grid(
  p_f05_grob_labeled, p_nnse_grob_labeled,  # Top row with labels
  island_heatmap_f05_grob_labeled, hist_f05a_grob_labeled,  # Bottom row with labels
  ncol = 2, nrow = 2,  # 2x2 grid layout
  rel_widths = c(1, 1),  # Equal widths for the columns
  rel_heights = c(1, 1)  # Equal heights for the rows
)


# Display the labeled 2x2 plot
print(fig2x2_labeled)

pdf(file   = "results/paper_figs/subset_heatmap_hist_v2.pdf",
    width  = 11,    # inches
    height = 9,
    family = "Helvetica"   # or another installed font
)
fig2x2_labeled
dev.off()

### ---- Fig. 4: ----

# Fig. 4a,b,c are jaccard_isl_f05
# Fig. 4d is cor_plot_dif_isl_f05

# helpers

# Remove any text/label annotation layers (e.g., annotate("text", ...))
drop_text_layers <- function(p) {
  if (!length(p$layers)) return(p)
  is_text_layer <- vapply(
    p$layers,
    function(L) any(grepl("GeomText|GeomLabel", class(L$geom), ignore.case = TRUE)),
    logical(1)
  )
  p$layers <- p$layers[!is_text_layer]
  p
}

# Compute "r = ..., p = ..." once
rp_text <- function(x, y, data, digits_r = 2, digits_p = 3) {
  ct <- cor.test(data[[y]], data[[x]], use = "complete.obs", method = "pearson")
  r_value <- round(unname(ct$estimate), digits_r)
  p_value <- if (ct$p.value < 1e-3) formatC(ct$p.value, format = "e", digits = 2)
  else formatC(ct$p.value, format = "f", digits = digits_p)
  paste0("r = ", r_value, ", p = ", p_value)
}

# Add a centered, plain header ABOVE the plot area and RESERVE space for it
# so headers don't overlap the plot and can't be clipped.
# Add a centered, plain header ABOVE the plot area and RESERVE space for it
# Works with older cowplot (no `padding` arg)
# Works with older cowplot: no padding=, no clip=
add_center_header <- function(p, header_text, size = 11, header_height = 0.12) {
  # Hide facet strips; we provide our own header
  p <- p + theme(strip.text = element_blank())
  
  # A tiny plot that only renders the centered header text
  header <- ggdraw() +
    draw_label(
      label = header_text,
      x = 0.5, y = 0.4,         # centered, near the top
      hjust = 0.5, vjust = 1,
      fontface = "plain", size = size
    )
  
  # Stack header above the plot; reserve vertical space via rel_heights
  plot_grid(
    header, p,
    ncol = 1,
    rel_heights = c(header_height, 1),
    align = "v"
  )
}

# plotting function (Jaccard)

make_facet_scatter_plot2 <- function(
    data,
    evaluator = "f05_score",
    pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
    names_to = "jaccard_type",
    values_to = "jaccard_value",
    facet_scales = "free_x",
    tme = theme_minimal()  # pass your theme here
) {
  df_long <- data %>%
    pivot_longer(cols = all_of(pivot_cols), names_to = names_to, values_to = values_to)
  
  facet_labels <- c(
    jaccard_edges        = "Interaction overlap",
    jaccard_plants       = "Plants overlap",
    jaccard_pollinators  = "Pollinators overlap"
  )
  
  ggplot(df_long, aes_string(x = values_to, y = evaluator)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "thistle") +
    facet_wrap(as.formula(paste("~", names_to)),
               scales = facet_scales,
               labeller = as_labeller(facet_labels)) +
    scale_x_continuous(labels = number_format(accuracy = 0.02)) +
    labs(x = "Jaccard similarity", y = "F0.5 score") +     # <-- (3) nice axis titles
    tme +                                            # your theme
    theme(
      strip.text   = element_text(size = 12, face = "plain"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks   = element_line(color = "black"),
      axis.text.x  = element_text(size = 14)
    )
}

# build panels

# (a) Pollinators
p_a_core <- make_facet_scatter_plot2(
  data = canary_results_jaccard,
  evaluator = "f05_score",
  pivot_cols = "jaccard_pollinators",
  facet_scales = "free_x",
  tme = tme
)
label_a <- paste0("Pollinator overlap: ", rp_text("jaccard_pollinators", "f05_score", canary_results_jaccard))
p_a <- add_center_header(p_a_core, label_a, size = 12)

# (b) Plants — keep y ticks; drop only the y-axis title (right column)
p_b_core <- make_facet_scatter_plot2(
  data = canary_results_jaccard,
  evaluator = "f05_score",
  pivot_cols = "jaccard_plants",
  facet_scales = "free_x",
  tme = tme
) + theme(axis.title.y = element_blank())
label_b <- paste0("Plant overlap: ", rp_text("jaccard_plants", "f05_score", canary_results_jaccard))
p_b <- add_center_header(p_b_core, label_b, size = 12)

# (c) Edges
p_c_core <- make_facet_scatter_plot2(
  data = canary_results_jaccard,
  evaluator = "f05_score",
  pivot_cols = "jaccard_edges",
  facet_scales = "free_x",
  tme = tme
)
label_c <- paste0("Interaction overlap: ", rp_text("jaccard_edges", "f05_score", canary_results_jaccard))
p_c <- add_center_header(p_c_core, label_c, size = 12)

# (d) Correlation — remove in-panel stats; header carries the stats
p_d_core <- cor_plot_dif_isl_f05 +
  tme +
  theme(
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 14)   # <— shrink x tick labels here
  )
p_d_core <- drop_text_layers(p_d_core)   # strip annotate("text", ...) if present
label_d <- paste0("Geographic distance: ", rp_text("distance_km", "f05_score", result_summary_island_dif))
p_d <- add_center_header(p_d_core, label_d, size = 12)

# arrange 2×2

isl_jaccard_distance <- plot_grid(
  p_a, p_b, p_c, p_d,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  ncol = 2,
  align = "hv",
  axis  = "tblr",
  label_size = 13,
  # raise the labels so headers sit below them (increase if your headers still touch)
  label_y = c(0.96, 0.96, 0.96, 0.96),
  label_x = c(0, 0, 0, 0),
  rel_widths = c(1,0.95,1,0.95)
  
)

isl_jaccard_distance

pdf(
  file   = "results/paper_figs/isl_jaccard_distance.pdf",
  width  = 8,    # inches
  height = 8,
  family = "Helvetica"   # or another installed font
)
print(isl_jaccard_distance)
dev.off()     # close the file

# ---- save results summary ----
# results_summary_islands <- summary(result_summary)
# 
# save_results_summary_from_summarytable(
#   results_summary_islands,
#   best_discrete_threshold = best_discrete_threshold,
#   out_csv = "results/predictions_island_scale_summary.csv"
# )
