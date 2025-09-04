# ---- Predicting interactions across space with SVD ----
# this pipeline allows us to predict missing links using the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# here we focus on island scale, but there is a section for comparison between scales.
# stages are according to the pipeline figure (Fig. 1).
### code for publication ###
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
library(softImpute)
library(ecodist)
library(rstatix)
library(ggrepel)

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- parameters ----
emln_id <- 60 # Canary Islands pollination system from Trøjelsgaard et al. 2015
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
  
  # Extract the reconstructed P matrix from C_reconstructed
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
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

# transforming raw predictions for binary evaluation
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# calculate balanced per-species f1 score
# 2) A helper that, given one species’ data.frame, 
#    does one draw of a balanced F₁
compute_balanced_f1 <- function(df_sp) {
  # split positives / negatives
  pos <- df_sp %>% filter(original_binary == 1)
  neg <- df_sp %>% filter(original_binary == 0)
  
  n_pos <- nrow(pos)
  n_neg <- nrow(neg)
  
  # if either class is missing, we can’t compute F1
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  
  # sample size = the smaller of the two
  n <- min(n_pos, n_neg)
  
  # draw one balanced sample
  samp_pos <- pos %>% sample_n(n)
  samp_neg <- neg %>% sample_n(n)
  samp     <- bind_rows(samp_pos, samp_neg)
  
  # compute TP, FP, FN
  TP <- sum(samp$original_binary == 1 & samp$predicted_bin_sigm == 1)
  FP <- sum(samp$original_binary == 0 & samp$predicted_bin_sigm == 1)
  FN <- sum(samp$original_binary == 1 & samp$predicted_bin_sigm == 0)
  
  # precision, recall
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  
  # F1
  if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    return(NA_real_)
  } else {
    return(2 * precision * recall / (precision + recall))
  }
}

# function to extract island names (removes "_site_X")
extract_island <- function(name) {
  gsub("_site_[12]", "", name)
}

# for MRM test
sym_average <- function(m) {
  mm <- m
  for(i in 1:nrow(mm)) for(j in 1:ncol(mm)) {
    if(i < j && !is.na(m[i,j]) && !is.na(m[j,i])) {
      avg       <- mean(c(m[i,j], m[j,i]))
      mm[i,j]   <- avg
      mm[j,i]   <- avg
    }
  }
  mm
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
      axis.text.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(values = custom_colors) + tme
}

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

plot_f1_nnse_vs_size_free_both <- function(data) {
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
        evaluator    = c(f1_score = "F1 score", nnse = "NNSE"),
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

plot_f1_nnse_vs_density_free_both <- function(data) {
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
        evaluator    = c(f1_score = "F1 score", nnse = "NNSE"),
        measure_type = c(density_P  = "Density of matrix P",
                         density_C  = "Density of matrix C")
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
      name   = "Network density",
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
      strip.text.x       = element_text(size = 12),
      strip.text.y.left  = element_text(size = 14, face = "bold", angle = 90),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks         = element_line(color = "black"),
      strip.background   = element_blank(),
      panel.spacing.x    = unit(0.7, "cm")
      
    )
}

make_facet_scatter_plot <- function(data,
                                    evaluator = "f1_score", 
                                    pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                    names_to = "jaccard_type", 
                                    values_to = "jaccard_value",
                                    x_lab = "Jaccard similarity",
                                    y_lab = "F1 score",
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
                                       evaluator = "f1_score",
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
                              y_axis_label = "F1 score",
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
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                            NA, 
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
      

      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
      
      sum(is.na(C))
      
      ### ---- b. prediction with SVD ----
      k_values <- 2
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
### ---- Fig. S7: selecting optimal threshold ----
# select the threshold for classifying links as 1s or 0s based on balance between f1 and balanced accuracy

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
               values_to = "value") %>%
  mutate(metric = recode(metric,
                         specificity       = "Specificity",
                         precision         = "Precision",
                         recall            = "Recall",
                         f1_score          = "F1 score",
                         balanced_accuracy = "Balanced accuracy",
                         mcc               = "MCC"
  ))

# 4) plot
optimal_threshold <- ggplot(df_avg, aes(threshold, value, color = metric)) +
  geom_line(size = 1) +
  labs(
    x     = "Probability threshold",
    y     = "Average metric",
    color = "Metric"
  ) +
  scale_color_brewer(palette = "Pastel2") +
  tme

# pdf(
#   file   = "optimal_threshold.pdf",
#   width  = 5,    # inches
#   height = 4,
#   family = "Helvetica"   # or another installed font
# )
# print(optimal_threshold)
# dev.off()     # close the file


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

df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

### ---- Fig. S8: predicted vs. observed weights ----
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
  theme_minimal(base_size = 14) + tme

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
    nse  = mean(nse,  na.rm = TRUE),
    nnse = mean(nnse, na.rm = TRUE)
  ) %>%
  ungroup()

head(result_summary)
summary(result_summary) # result_summary includes evaluation results across all iterations for each combination of islands 

### ---- Fig. 2b: distribution of evaluators with/without external data ----
# this analysis shows us if predictions made using added information from other locations (off-diagonals in layer-to-layer predictions, as a heatmap) is any better than not adding any information (cases on the diagonal)

result_summary <- result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Single location",
    train_layer != test_layer ~ "Added location"
  ))

custom_colors <- c("Single location" = "steelblue",
                   "Added location" = "thistle")

# plot the histogram: 

hist_f1a <- plot_hist(result_summary, metric = "f1_score", 
                     y_axis_label = "Count of instances",
                     x_axis_label = "F1 score") + 
  scale_y_continuous(labels = scales::number_format(accuracy = 1.0))

hist_f1a # Fig. 2b

# # Base‐R PDF device
# pdf(
#   file   = "hist_f1a.pdf",
#   width  = 5,    # inches
#   height = 4,
#   family = "Helvetica"   # or another installed font
# )
# print(hist_f1a)
# dev.off()     # close the file

# stats
# run t-test via formula interface
t_test_f1 <- t.test(f1_score ~ layer_comparison, 
                    data       = result_summary,
                    var.equal  = FALSE)  # Welch’s test

# 3. Print the full test
print(t_test_f1)

# do we need welch/wilcoxon?
# first normality check
result_summary %>%
  group_by(layer_comparison) %>%
  shapiro_test(f1_score)

# variance check
result_summary %>% levene_test(f1_score ~ layer_comparison)
# all is good, we can use t-test.

# 4. Extract just the numbers you want
t_stat <- unname(t_test_f1$statistic)
df_val <- unname(t_test_f1$parameter)
p_val  <- t_test_f1$p.value

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
    C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                          NA, 
                                          C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
    
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

# summerize (table ST1)
# first add layer names
# add distances to the main table
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# create a new grouped tibble
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

#### ---- Fig. 5: correlate network size with evaluators ----

df_netsize <- result_summary %>%
  select(f1_score, nnse, size_P, density_P, size_C, density_C) %>%
  pivot_longer(
    cols = c(size_P, density_P, size_C, density_C),
    names_to = "measure_type",
    values_to = "measure_value"
  )


df_f1_nnse_size <- result_summary %>%
  select(f1_score, nnse, size_P, size_C) %>%
  pivot_longer(cols = c(size_P, size_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f1_score, nnse), names_to = "evaluator", values_to = "evaluator_value")

netsize_f1_nnse <- plot_f1_nnse_vs_size_free_both(df_f1_nnse_size) + tme # Fig. 5
netsize_f1_nnse

# pdf(
#   file   = "netsize_f1_nnse.pdf",
#   width  = 6,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(netsize_f1_nnse)
# dev.off()     # close the file

#### ---- Fig. S12: density ----

df_f1_nnse_density <- result_summary %>%
  select(f1_score, nnse, density_P, density_C) %>%
  pivot_longer(cols = c(density_P, density_C), names_to = "measure_type", values_to = "measure_value") %>%
  pivot_longer(cols = c(f1_score, nnse), names_to = "evaluator", values_to = "evaluator_value")

netdensity_f1_nnse <- plot_f1_nnse_vs_density_free_both(df_f1_nnse_density) + tme
netdensity_f1_nnse

# pdf(
#   file   = "netdensity_f1_nnse.pdf",
#   width  = 6,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(netdensity_f1_nnse)
# dev.off()     # close the file

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

#### ---- Fig. 4: plot Jaccard ----
# we filter only pairs of different islands for this analysis, since Jaccard index for the same island is 1 for all islands
canary_results_jaccard <- result_summary %>%
  filter(train_layer != test_layer)

jaccard_isl_f1 <- make_facet_scatter_plot(data = canary_results_jaccard, 
                                        evaluator = "f1_score",
                                        pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                        x_lab = "Jaccard similarity",
                                        y_lab = "F1 score",
                                        facet_scales = "free_x")
jaccard_isl_f1 # Fig. 4

# # Base‐R PDF device
# pdf(
#   file   = "jaccard_isl_f1.pdf",
#   width  = 7,    # inches
#   height = 3.5,
#   family = "Helvetica"   # or another installed font
# )
# print(jaccard_isl_f1)
# dev.off()     # close the file

jaccard_isl_nnse <- make_facet_scatter_plot(data = canary_results_jaccard, 
                                           evaluator = "nnse",
                                           pivot_cols = c("jaccard_pollinators", "jaccard_plants", "jaccard_edges"),
                                           x_lab = "Jaccard similarity",
                                           y_lab = "NNSE",
                                           facet_scales = "free_x")
jaccard_isl_nnse

### ---- partner fidelity correlation with evaluators ----

# filter existing interactions in cases where train = test layer (we want to calculate the average similarity in partner composition of each species in different locations)
df_fidelity <- df %>% filter(train_layer == test_layer) %>% 
  filter(original_links != 0) %>% filter(itr == 1) # all iterations are the same for existing links that were not removed

#### ---- plants fidelity ----
# for each plant (node_from) and layer, gather the pollinators (node_to).
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
      # If the plant is only in one layer, there are no pairs, so return NA. select only species the occur in at least 3 islands
      if (n_layers < 3) {
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

# we filter out single-island predictions since what is relevant here is prediction across space
df_off <- df %>% 
  filter(train_layer != test_layer) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

df_merged <- df_off %>%
  left_join(df_sorensen, by = "node_from")

#### ---- pollinators fidelity ----
df_pollinators <- df_fidelity %>%
  distinct(node_to, train_layer, node_from) %>%
  group_by(node_to, train_layer) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# 2. For each pollinator, compute mean Sorensen similarity across all pairs of layers.
df_sorensen_pollinators <- df_pollinators %>%
  group_by(node_to) %>%
  summarise(
    mean_sorensen_pollinators = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 3) {
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
  left_join(df_sorensen_pollinators, by = "node_to") # joint plant and pollinator fidelities

#### ---- per-species f1 ----
# how many removed 1s and 0s do we have for each species?
# For plants:
removed_count_plants <- df_merged %>%
  filter(removed == 1) %>% 
  group_by(plant = node_from) %>%
  summarise(
    n_removed_0    = sum(original_binary == 0, na.rm = TRUE),
    n_removed_1    = sum(original_binary == 1, na.rm = TRUE),
  ) %>%
  arrange(plant)

# For pollinators:
removed_count_polls <- df_merged %>%
  filter(removed == 1) %>% 
  group_by(pollinator = node_to) %>%
  summarise(
    n_removed_0    = sum(original_binary == 0, na.rm = TRUE),
    n_removed_1    = sum(original_binary == 1, na.rm = TRUE),
  ) %>%
  arrange(pollinator)

# calculate balanced f1 for each species by sampling its removed links 50 times in a balanced subset and calculating f1 for each
df_clean <- df_merged %>%
  filter(removed == 1) %>% 
  filter(!is.na(original_binary), !is.na(predicted_bin_sigm))

set.seed(42)

# 3) For pollinators (node_to), bootstrap 50× per species
df_balanced_f1_polls <- df_clean %>%
  group_by(species = node_to) %>%
  group_modify(~ {
    f1s <- replicate(50, compute_balanced_f1(.x))
    tibble(
      mean_f1 = mean(f1s, na.rm = TRUE),
      sd_f1   = sd(  f1s, na.rm = TRUE)
    )
  }) %>%
  ungroup()

# 4) Do the same for plants (node_from)
df_balanced_f1_plants <- df_clean %>%
  group_by(species = node_from) %>%
  group_modify(~ {
    f1s <- replicate(50, compute_balanced_f1(.x))
    tibble(
      mean_f1 = mean(f1s, na.rm = TRUE),
      sd_f1   = sd(  f1s, na.rm = TRUE)
    )
  }) %>%
  ungroup()

# join them together
df_for_plot_balanced <- bind_rows(
  df_balanced_f1_polls %>%
    left_join(df_sorensen_pollinators, by = c("species" = "node_to")) %>%
    mutate(avg_sorensen_plants = NA_real_),
  
  df_balanced_f1_plants %>%
    left_join(df_sorensen, by = c("species" = "node_from")) %>%
    mutate(avg_sorensen_pollinators = NA_real_)
)

#### ---- Fig. S13: plot fidelity ----
balanced_f1_fidelity <- make_full_correlation_plot(
  data            = df_for_plot_balanced,
  evaluator       = "mean_f1",
  pollinator_x    = "mean_sorensen_pollinators",
  plant_x         = "mean_sorensen_plants",
  shared_x_lab    = "Mean partner fidelity (Sørensen)",
  shared_y_lab    = "F1 score per species"
)

balanced_f1_fidelity

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
  summarise(count_never_observed = n(), .groups = "drop")

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
    y = "Number of predicted, \nnon-observed interactions",
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

final_plot <- combine_plots(plant_degree, poll_degree) # fix error
grid::grid.newpage()
grid::grid.draw(final_plot) # Fig. 6c

# pdf("degree_unobserved_links.pdf", width = 10, height = 7)
# grid::grid.draw(final_plot)
# dev.off()

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

# reset the levels for the species factors
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

### ---- Fig. S10: difference in links predicted with/without external data ----
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

map_missing_links_diags_offs <- ggplot(df_plot, aes(x = node_to, y = node_from)) +
  
  # allow a second fill scale
  new_scale_fill() +
    geom_tile(
    data  = filter(df_plot, !is.na(sigm_cat)),
    aes(fill = sigm_cat),
    alpha = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "offs↑ only" = "salmon",
      "diag↑ only" = "plum3",
      "both↑"       = "lightsteelblue"
    ),
    na.value = NA,
    name   = "Difference in \nprediction approach",
    labels = c(
      "offs↑ only" = "Predicted only by \nadding external location",
      "diag↑ only" = "Predicted only by \nsingle location",
      "both↑"       = "Predicted by both approaches"
    )
  ) +
  
  # tidy up
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 8),
    legend.position = "bottom",
    legend.box      = "vertical"
  ) +
  
  # your italic‐species labels and extra theme element
  scale_y_discrete(
    labels = function(x) lapply(strsplit(x, "_"), function(y) {
      bquote(italic(.(paste(y, collapse = " "))))
    })
  ) +
  tme

map_missing_links_diags_offs

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
### ---- Fig. 3b: pie chart ----
# 3. Make the pie
pie_chart <- ggplot(df_counts, aes(x = "", y = n, fill = sigm_cat)) +
  geom_col(width = 1, color = "white") +      # white border between slices
  coord_polar(theta = "y") +                  # convert bar → pie
  scale_fill_manual(values = my_cols,
                    labels = new_labels) +
  theme_void() +                              # remove axes/background
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 17),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position  = "bottom",
    legend.direction = "vertical"
  ) +
  #labs(title = "Interactions by Significance Category") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 7
  )

pie_chart
# pdf(
#   file   = "pie_chart.pdf",
#   width  = 7,    # inches
#   height = 7,
#   family = "Helvetica"   # or another installed font
# )
# print(pie_chart)
# dev.off()     # close the file

### ---- distance decay ----
#### ---- add distances and location names ----
distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL)

# proceed for both island scale and site scale and compare the trends

#### ---- island scale ----
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

#### ---- site scale ----
# add names and distances to the main table
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

# load site scale data
#result_site <- read_csv('canary_weighted_scaled_site_net_60_50_itr.csv') # from file
# or run the analysis at site scale
set.seed(42)

# Load matrices
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
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
result_site <- data.frame()

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
      # Ensure that existing entries are not overwritten; sum overlapping entries
      C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]),
                                            NA,
                                            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])


      # Apply biScale to center matrices
      C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)

      sum(is.na(C))

      ### ---- b. prediction with SVD ----
      k_values <- 2
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

    result_site <- rbind(
      result_site,
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

# # result_site includes predictions for all combinations of sites, 50 iterations of links withholding and prediction for each combination
# convert negatives to zeros
result_site <- result_site %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# create the evaluation table
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

#### ---- distance correlation with evaluators ----

# remove sites form within the same island - only use information from different islands for distance decay
result_summary_island_dif <- result_summary_island %>% filter(train_layer != test_layer)
result_summary_site_dif <- result_summary_site %>% filter(train_layer != test_layer)

# extract island names
result_summary_site_dif$train_island <- sub("^(\\w+).*", "\\1", result_summary_site_dif$train_layer_name)
result_summary_site_dif$test_island <- sub("^(\\w+).*", "\\1", result_summary_site_dif$test_layer_name)

# filter rows where island names are different
filtered_results <- result_summary_site_dif[result_summary_site_dif$train_island != result_summary_site_dif$test_island, ]

# plot
cor_plot_site_dif_f1 <- make_cor_plot(filtered_results, evaluator = "f1_score", extra_theme = tme)
cor_plot_dif_isl_f1  <- make_cor_plot(result_summary_island_dif, evaluator = "f1_score", extra_theme = tme)
distance_dif_plot_f1 <- combine_two_plots(cor_plot_site_dif_f1, cor_plot_dif_isl_f1)

# run the previous make_cor_plot again
# pdf(
#   file   = "distance_plot_f1_different_isl.pdf",
#   width  = 7,
#   height = 4,
#   family = "Helvetica"
# )
# 
# grid::grid.draw(distance_dif_plot_f1)
# 
# dev.off()

#### ---- MRM test: island scale ----

# build square matrices of f1 and distance
#    (layers must be in the same order for rows & cols)

# a) get a list of all unique layers
layers <- sort(unique(c(result_summary_island_dif$train_layer, result_summary_island_dif$test_layer)))

# b) initialize empty matrices
f1_mat      <- matrix(NA, nrow=length(layers), ncol=length(layers),
                      dimnames=list(layers, layers))
dist_mat_km <- f1_mat

# c) fill in each cell [i,j] with the corresponding f1_score and distance_km
for(i in layers) for(j in layers) {
  # subset rows where train=i and test=j
  sub <- result_summary_island_dif[result_summary_island_dif$train_layer==i & result_summary_island_dif$test_layer==j, ]
  if(nrow(sub)==1) {
    f1_mat[i,j]      <- sub$f1_score
    dist_mat_km[i,j] <- sub$distance_km
  }
}

# d) because MRM uses symmetric distance matrices, average [i,j] & [j,i]
f1_sym      <- sym_average(f1_mat)
dist_sym_km <- sym_average(dist_mat_km)

# e) convert to “dist” objects (lower triangle)
dist_f1      <- as.dist(f1_sym)
dist_km      <- as.dist(dist_sym_km)

# 4. run the MRM
#    — this will regress the F1‐distance matrix on the geographic–distance matrix
set.seed(42)   # for reproducibility of permutations
mrm_out <- MRM(dist_f1 ~ dist_km, nperm=999)

# 5. results
print(mrm_out)

#### ---- MRM for site scale ----
layers_site <- sort(unique(c(result_summary_site_dif$train_layer, result_summary_site_dif$test_layer)))

# initialize empty matrices
f1_mat_site      <- matrix(NA, nrow=length(layers_site), ncol=length(layers_site),
                      dimnames=list(layers_site, layers_site))
dist_mat_km_site <- f1_mat_site

# fill in each cell [i,j] with the corresponding f1_score and distance_km
for(i in layers_site) for(j in layers_site) {
  # subset rows where train=i and test=j
  sub <- result_summary_site_dif[result_summary_site_dif$train_layer==i & result_summary_site_dif$test_layer==j, ]
  if(nrow(sub)==1) {
    f1_mat_site[i,j]      <- sub$f1_score
    dist_mat_km_site[i,j] <- sub$distance_km
  }
}

# because MRM uses symmetric distance matrices, average [i,j] & [j,i]

f1_sym_site      <- sym_average(f1_mat_site)
dist_sym_km_site <- sym_average(dist_mat_km_site)

# convert to “dist” objects (lower triangle)
dist_f1_site      <- as.dist(f1_sym_site)
dist_km_site     <- as.dist(dist_sym_km_site)

# run the MRM
#    — this will regress the F1‐distance matrix on the geographic–distance matrix
set.seed(42)   # for reproducibility of permutations
mrm_out_site <- MRM(dist_f1_site ~ dist_km_site, nperm=999)

# results
print(mrm_out_site)

### ---- Fig. 2a heatmap ----

island_heatmap_f1 <- 
  ggplot(result_summary_island, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_island[result_summary_island$train_layer == result_summary_island$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added location", y = "Predicted location", fill = "F1 score") +
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

# pdf(
#   file   = "island_heatmap_f1.pdf",
#   width  = 6,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(island_heatmap_f1)
# dev.off()     # close the file

### ---- Fig. S1: compare scales ----

df_long <- bind_rows(
  result_summary_site   %>% mutate(scale = "Site"),
  result_summary_island %>% mutate(scale = "Island")
) %>%
  pivot_longer(
    cols      = c("f1_score", "nnse"),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = factor(metric, levels = c(
    "f1_score", "nnse"
  )))

# test for assumptions
# normality

df_long %>%
  group_by(metric, scale) %>%
  shapiro_test(value)   # Shapiro-Wilk normality test

# df_long %>%
#   ggqqplot(x = "value", facet.by = c("metric", "scale"))

# mostly, data is not normally distributed, so better use wilcoxon 

# plot
# pretty facet titles with units:
metric_labels <- c(
  f1_score          = "F1 score",
  nnse              = "NNSE"
)

nnse_f1_scales <- ggplot(df_long, aes(x = scale, y = value, fill = scale)) +
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
    switch        = "y"             # move the strip to the left side
  ) +
  stat_compare_means(
    method         = "wilcox.test",
    label          = "p.format",    # print the full p‐value
    p.format.args  = list(
      digits     = 2,               # two digits after decimal
      scientific = TRUE             # use e-notation for small p’s
    ),
    label.y        = Inf,
    vjust          = 1.5,
    label.x        = 1.45,
    tip.length     = 0.01,
    size           = 3.5              # adjust this for font size
  ) +
  scale_fill_manual(values = c("Site"   = "lightsteelblue2",
                               "Island" = "wheat2")) +
  labs(
    x = NULL,
    y = NULL                       # we’ll rely on the left‐side strips as “y‐titles”
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement       = "outside",           # draw strips outside the plot panel
    strip.text.y.left     = element_text(
      angle = 90,          # horizontal text
      face  = "bold",
      size  = 12
    ),
    axis.text.x           = element_blank(),
    axis.ticks.x          = element_blank(),
    legend.position       = "bottom"
  ) + tme

nnse_f1_scales

#Base‐R PDF device
# pdf(
#   file   = "scales_fig2.pdf",
#   width  = 7,    # inches
#   height = 4,
#   family = "Helvetica"   # or another installed font
# )
# print(nnse_f1_scales)
# dev.off()     # close the file
