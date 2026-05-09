
# libraries -------------
library(emln)
library(ggplot2)
library(softImpute)

# parameters -------------
emln_id <- 60 # Canary Islands pollination system from Trøjelsgaard et al. 2015
prop_ones_to_remove <- 0.2 # proportion of existing links to withhold
default_threshold <- 0.6 # threshold optimal for F0.5 score, calc in main script

# themes -----------
tme <-  theme(axis.text = element_text(size = 18, color = "black"),
              axis.title = element_text(size = 18, face = "bold"),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())


# functions -------------
# transforming raw predictions for binary evaluation
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

build_interaction_matrix <- function(data, layers_to_filter) {
  
  # unify layer format
  layers <- layers_to_filter # default, they are character
  if (is.numeric(layers_to_filter)) layers <- paste0("layer_", layers_to_filter)
  
  # Step 1: Filter rows based on specified layers
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
implement_impute <- function(C, k, lambda, P, remove_indices, zeros_to_remove_indices, P_original, back_trans_values = NULL) {
  if (!is.null(back_trans_values)) { 
    row_centers <- back_trans_values$row_centers
    col_centers <- back_trans_values$col_centers
  } # if it is NULL, then the values are already in global variables. not recommended, though.
  
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
    if (is.null(dim(zeros_to_remove_indices))) { # handles when remove_indices has only one row
      test_indices <- rbind(
        data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
        data.frame(row = zeros_to_remove_indices["row"], col = zeros_to_remove_indices["col"], label = 0)
      )
    } else{
      test_indices <- rbind(
        data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
        data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"], label = rep(0, nrow(zeros_to_remove_indices)))
      )
    }
    

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


plot_f05_nnse_vs_size_free_both <- function(data) {
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
  
  # panel labels
  panel_labels <- tibble::tribble(
    ~evaluator,  ~measure_type, ~panel_label,
    "f05_score", "size_C",   "(a)",
    "f05_score", "size_P",   "(b)",
    "nnse",      "size_C",   "(c)",
    "nnse",      "size_P",   "(d)"
  )
  
  ggplot(data, aes(x = measure_value, y = evaluator_value)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "salmon") +
    
    facet_grid(
      rows   = vars(evaluator),
      cols   = vars(measure_type),
      scales = "free",     # ← free both x and y per facet
      labeller = labeller(
        evaluator    = c(f05_score = "F0.5 score", nnse = "NNSE"),
        measure_type = c(size_P  = "Size of matrix P",
                         size_C  = "Size of matrix C")
      ),
      switch = "y"
    ) +
    
    geom_text(
      data        = cor_table,
      aes(label   = label_text),
      x           = Inf,
      y           = Inf,
      hjust       = 1.1,
      vjust       = 1.2,
      size        = 3.2,
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data        = panel_labels,
      aes(label   = panel_label),
      x           = -Inf,
      y           = Inf,
      hjust       = -0.4,
      vjust       = 1.4,
      size        = 5,
      fontface    = "bold",
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

plot_f05_nnse_vs_density_free_both <- function(data) {
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
  
  # panel labels
  panel_labels <- tibble::tribble(
    ~evaluator,  ~measure_type, ~panel_label,
    "f05_score", "density_C",   "(a)",
    "f05_score", "density_P",   "(b)",
    "nnse",      "density_C",   "(c)",
    "nnse",      "density_P",   "(d)"
  )
  
  ggplot(data, aes(x = measure_value, y = evaluator_value)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "salmon") +
    
    facet_grid(
      rows   = vars(evaluator),
      cols   = vars(measure_type),
      scales = "free",
      labeller = labeller(
        evaluator    = c(f05_score = "F0.5 score", nnse = "NNSE"),
        measure_type = c(density_P  = "Connectance of matrix P",
                         density_C  = "Connectance of matrix C")
      ),
      switch = "y"
    ) +
    
    geom_text(
      data        = cor_table,
      aes(label   = label_text),
      x           = Inf,
      y           = Inf,
      hjust       = 1.1,
      vjust       = 1.2,
      size        = 3.2,
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data        = panel_labels,
      aes(label   = panel_label),
      x           = -Inf,
      y           = Inf,
      hjust       = -0.4,
      vjust       = 1.4,
      size        = 5,
      fontface    = "bold",
      inherit.aes = FALSE
    ) +
    
    scale_x_continuous(
      name   = "Network connectance",
      breaks = scales::breaks_width(0.02),
      labels = scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    
    scale_y_continuous(
      name   = NULL,
      expand = expansion(mult = c(0.05, 0.1))
    ) +
    
    theme_minimal() +
    theme(
      strip.placement    = "outside",
      strip.text.x       = element_text(size = 12),
      strip.text.y.left  = element_text(size = 14, face = "bold", angle = 90),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.ticks         = element_line(color = "black"),
      strip.background   = element_blank(),
      panel.spacing.x    = unit(0.7, "cm"),
      axis.text.x        = element_text(size = 10),
      axis.text.y        = element_text(size = 10)
    )
}


load_and_mold_data_for_prediction <- function(data_id){
  d <- load_emln(data_id)
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
  
  # set new layer names using the old ones
  aggregated_df <- aggregated_df %>% 
    separate_wider_delim(layer_from, delim = "_", names = c("t", "l1", "l2"), cols_remove = FALSE) %>%
    mutate(island_id = paste0("layer_", as.numeric(l2)/2))  %>%
    mutate(layer_from = island_id, layer_to = island_id)%>%
    select(layer_from, node_from, layer_to, node_to, weight, type)
  
  return(aggregated_df)
}

get_island_names_with_layer_indexes <- function() {
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
  
}

# save results summary
save_results_summary_from_summarytable <- function(summary_table_obj,
                                                   best_discrete_threshold,
                                                   out_csv) {
  # Convert the table (chr matrix) to a character matrix
  mat <- as.matrix(summary_table_obj)
  
  # Column names are your variables (may have leading spaces)
  vars <- trimws(colnames(mat))
  if (is.null(vars) || length(vars) == 0) {
    stop("No column names found on the summary table object.")
  }
  
  # Helper: extract the first numeric value from strings like "Min.   :60"
  extract_num <- function(x) {
    x <- as.character(x)
    m <- regexpr("-?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?", x)
    ifelse(m == -1, NA_real_, as.numeric(regmatches(x, m)))
  }
  
  # Your table rows are typically: Min, 1st Qu, Median, Mean, 3rd Qu, Max
  row_labels <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
  
  # Parse values column-by-column
  vals <- lapply(seq_along(vars), function(j) extract_num(mat[, j]))
  vals <- do.call(rbind, vals)  # variables x 6
  colnames(vals) <- row_labels
  
  df_out <- data.frame(
    Variable = vars,
    vals,
    best_discrete_threshold = best_discrete_threshold,
    row.names = NULL,
    check.names = FALSE
  )
  
  # Ensure output directory exists
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  
  # Save
  write.csv(df_out, out_csv, row.names = FALSE)
  message("Saved: ", out_csv)
  
  invisible(df_out)
}