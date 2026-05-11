# ---- Temporal host-parasite interaction prediction ----
# this pipeline allows us to predict missing links using the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# here we focus on year scale (temporal scale), comparing different time periods.
# stages are according to the pipeline figure (Fig. 1).
## ---- Load libraries ----
library(tidyverse)
library(scales)
library(cowplot)
library(ggnewscale)
library(softImpute)
library(pROC)
library(PRROC)
library(readxl)

source("code/common.R")

## ---- Parameters ----
n_sim <- 50
set.seed(42)

## ---- Functions ----

# Override implement_impute: temporal prediction loop uses global variables
# (row_centers, col_centers, P, remove_indices, zeros_to_remove_indices, P_original)
implement_impute <- function(C, k, lambda) {
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  C_reconstructed <- softImpute::complete(C, fit)
  C_reconstructed_orig <- C_reconstructed +
    outer(row_centers[rownames(C)], col_centers[colnames(C)], "+")
  P_reconstructed <- C_reconstructed_orig[rownames(P), colnames(P)]

  if (is.null(dim(remove_indices))) {
    test_indices <- rbind(
      data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"],
                 label = rep(0, nrow(zeros_to_remove_indices)))
    )
  } else {
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"],
                 label = rep(1, nrow(remove_indices))),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"],
                 label = rep(0, nrow(zeros_to_remove_indices)))
    )
  }

  test_rows      <- rownames(P)[test_indices$row]
  test_cols      <- colnames(P)[test_indices$col]
  original_links <- P_original[cbind(test_rows, test_cols)]
  predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]

  results <- data.frame(k = k, lambda = lambda,
                        original_links   = original_links,
                        predicted_values = predicted_values,
                        node_to  = test_rows,
                        node_from = test_cols,
                        removed  = 1)

  all_edges <- expand.grid(
    node_to = rownames(P), node_from = colnames(P),
    k = k, lambda = lambda,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  all_edges$original_links <- mapply(
    function(r, c) P_original[r, c],
    all_edges$node_to, all_edges$node_from
  )
  removed_edges_idx <- data.frame(node_to = test_rows, node_from = test_cols,
                                  stringsAsFactors = FALSE)
  not_removed <- all_edges[
    !paste(all_edges$node_to, all_edges$node_from) %in%
      paste(removed_edges_idx$node_to, removed_edges_idx$node_from), ]
  not_removed$removed          <- 0
  not_removed$predicted_values <- NA
  not_removed$k      <- k
  not_removed$lambda <- lambda

  list(results = results, not_removed = not_removed)
}

# Scatter plot of F0.5 vs. Jaccard (used in Fig. S22 panels a, b, c)
make_facet_scatter_plot2 <- function(
    data,
    evaluator    = "f05_score",
    pivot_cols   = c("jaccard_parasites", "jaccard_hosts", "jaccard_edges"),
    names_to     = "jaccard_type",
    values_to    = "jaccard_value",
    facet_scales = "free_x",
    tme          = theme_minimal()) {

  df_long <- data %>%
    pivot_longer(cols = all_of(pivot_cols), names_to = names_to, values_to = values_to)

  facet_labels <- c(
    jaccard_edges       = "Interaction overlap",
    jaccard_hosts       = "Hosts overlap",
    jaccard_parasites   = "Parasites overlap",
    jaccard_plants      = "Hosts overlap",
    jaccard_pollinators = "Parasites overlap"
  )

  ggplot(df_long, aes_string(x = values_to, y = evaluator)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "thistle") +
    facet_wrap(as.formula(paste("~", names_to)),
               scales   = facet_scales,
               labeller = as_labeller(facet_labels)) +
    scale_x_continuous(labels = number_format(accuracy = 0.02)) +
    labs(x = "Jaccard similarity", y = "F0.5 score") +
    tme +
    theme(
      strip.text   = element_text(size = 12, face = "plain"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks   = element_line(color = "black"),
      axis.text.x  = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
}

# Correlation scatter plot (used in Fig. S22 panel d)
make_cor_plot <- function(data, evaluator,
                          distance_col = "time_difference",
                          x_lab        = "Time difference (years)",
                          y_lab        = NULL,
                          extra_theme  = NULL) {
  if (is.null(y_lab)) y_lab <- evaluator
  correlation <- cor.test(data[[evaluator]], data[[distance_col]],
                          use = "complete.obs", method = "pearson")
  r_value <- round(correlation$estimate, 2)
  p_value <- ifelse(
    correlation$p.value < 0.001,
    formatC(correlation$p.value, format = "e", digits = 2),
    formatC(correlation$p.value, format = "f", digits = 3)
  )
  label_text <- paste0("r = ", r_value, ", p = ", p_value)
  plot <- ggplot(data, aes_string(x = distance_col, y = evaluator)) +
    geom_point(color = "salmon2", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +
    labs(x = x_lab, y = y_lab) +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.1, size = 3.5, color = "black")
  if (!is.null(extra_theme)) plot <- plot + extra_theme
  return(plot)
}

# Helper: remove annotate("text",...) layers (used in Fig. S22 panel d assembly)
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

# Helper: compute "r = ..., p = ..." string for panel headers in Fig. S22
rp_text <- function(x, y, data, digits_r = 2, digits_p = 3) {
  ct      <- cor.test(data[[y]], data[[x]], use = "complete.obs", method = "pearson")
  r_value <- round(unname(ct$estimate), digits_r)
  p_value <- if (ct$p.value < 1e-3) formatC(ct$p.value, format = "e", digits = 2)
             else formatC(ct$p.value, format = "f", digits = digits_p)
  paste0("r = ", r_value, ", p = ", p_value)
}

# Helper: add a centered text header above a panel in Fig. S22
add_center_header <- function(p, header_text, size = 11, header_height = 0.12) {
  p <- p + theme(strip.text = element_blank())
  header <- ggdraw() +
    draw_label(label = header_text, x = 0.5, y = 0.4,
               hjust = 0.5, vjust = 1, fontface = "plain", size = size)
  plot_grid(header, p, ncol = 1, rel_heights = c(header_height, 1), align = "v")
}


## ---- 1. Load and preprocess data ----
hp_network <- readxl::read_excel("data/41559_2017_BFs415590170101_MOESM36_ESM.xlsx")

# Aggregate individual hosts to species level
hp_network <- hp_network %>%
  group_by(YearCollected, Host) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Reshape to long format matching the standard aggregated_df structure
aggregated_df <- hp_network %>%
  mutate(layer_from = paste0("layer_", YearCollected),
         layer_to   = paste0("layer_", YearCollected)) %>%
  pivot_longer(
    cols      = -c(YearCollected, Host, layer_from, layer_to),
    names_to  = "node_to",
    values_to = "weight"
  ) %>%
  rename(node_from = Host) %>%
  filter(weight > 0) %>%
  mutate(type = "host-parasite") %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# Rename layers to sequential numbers (layer_1, layer_2, ...)
unique_years  <- sort(unique(hp_network$YearCollected))
year_mapping  <- setNames(paste0("layer_", seq_along(unique_years)),
                          paste0("layer_", unique_years))
aggregated_df <- aggregated_df %>%
  mutate(layer_from = year_mapping[layer_from],
         layer_to   = year_mapping[layer_to])

# Reverse mapping: layer_1 -> 1982, etc.
layer_to_year <- setNames(unique_years, paste0("layer_", seq_along(unique_years)))

# Data frame for joining layer numbers to year labels
layer_year_df <- data.frame(
  layer_num = seq_along(layer_to_year),
  year      = as.character(layer_to_year),
  stringsAsFactors = FALSE
)

num_layers   <- length(unique(aggregated_df$layer_from))
results_file <- "results/hp_analysis/predictions_year_scale.rds"


## ---- 2. Load or compute predictions ----
combined_results <- data.frame()

if (file.exists(results_file)) {
  print("Existing results file found — reading the file and proceeding to analysis")
  combined_results <- readRDS(results_file)
  print("finished loading prediction results")

} else {
  for (layers_to_train in 1:num_layers) {
    for (layer_to_predict in 1:num_layers) {
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))

      A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
      P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)

      num_1_to_remove <- floor(sum(P > 0, na.rm = TRUE) * prop_ones_to_remove)
      ones_in_P       <- which(P > 0, arr.ind = TRUE)
      num_0_to_remove <- num_1_to_remove
      prop_0_removed  <- num_0_to_remove / sum(P == 0, na.rm = TRUE)
      zeros_in_P      <- which(P == 0, arr.ind = TRUE)

      bootstrapping_results <- NULL
      P_original <- P

      for (i in 1:n_sim) {
        remove_indices          <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
        P[remove_indices]       <- NA
        zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
        P[zeros_to_remove_indices] <- NA

        all_row_ids <- unique(c(rownames(A), rownames(P)))
        all_col_ids <- unique(c(colnames(A), colnames(P)))
        C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                    dimnames = list(all_row_ids, all_col_ids))
        C[rownames(A), colnames(A)] <- A

        if (layers_to_train != layer_to_predict) {
          C[rownames(P), colnames(P)] <- ifelse(
            is.na(C[rownames(P), colnames(P)]), NA,
            C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)]
          )
        } else {
          C[rownames(P), colnames(P)] <- ifelse(
            is.na(C[rownames(P), colnames(P)]), NA,
            (C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)]) / 2
          )
        }

        C <- biScale(C, row.center = TRUE, col.center = TRUE,
                     row.scale = FALSE, col.scale = FALSE)
        row_centers <- attr(C, "biScale:row")$center
        col_centers <- attr(C, "biScale:column")$center

        k_values      <- c(2, 5, 10)
        lam0          <- lambda0(C)
        lambda_values <- c(1, 5, 50, 100, lam0)

        results <- data.frame(k = integer(), lambda = numeric(),
                              original_links = numeric(), predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL

        for (k in k_values) {
          for (lambda in lambda_values) {
            r <- implement_impute(C, k, lambda)
            r$results$input_lambda     <- lambda
            r$not_removed$input_lambda <- lambda
            results         <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }

        complete_edges_all     <- rbind(results, not_removed_all)
        complete_edges_all$itr <- i
        bootstrapping_results  <- rbind(bootstrapping_results, complete_edges_all)
        P <- P_original
      }

      combined_results <- rbind(
        combined_results,
        cbind(
          data.frame(
            train_layer         = layers_to_train,
            test_layer          = layer_to_predict,
            prop_ones_removed   = prop_ones_to_remove,
            amount_of_removed_1 = num_1_to_remove,
            amount_of_removed_0 = num_0_to_remove,
            prop_0_removed      = prop_0_removed
          ),
          bootstrapping_results
        )
      )
    }
  }
  saveRDS(combined_results, file = results_file)
}


## ---- 3. Filter predictions ----
combined_results <- combined_results %>%
  filter(k == 2) %>%
  filter(!(input_lambda %in% c(1, 5, 50, 100)))

df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))


## ---- 4. Threshold selection ----
thresholds <- seq(0, 1, by = 0.1)

df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1, 0)
  )

df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%
  mutate(predicted_bin = if_else(predicted_prob > threshold, 1, 0)) %>%
  group_by(train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity       = TN / (TN + FP),
    precision         = TP / (TP + FP),
    recall            = TP / (TP + FN),
    f05_score         = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc  = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)
  ) %>%
  ungroup() %>%
  group_by(train_layer, test_layer, threshold) %>%
  summarise(
    TP                = mean(TP,                na.rm = TRUE),
    FN                = mean(FN,                na.rm = TRUE),
    TN                = mean(TN,                na.rm = TRUE),
    FP                = mean(FP,                na.rm = TRUE),
    specificity       = mean(specificity,       na.rm = TRUE),
    precision         = mean(precision,         na.rm = TRUE),
    recall            = mean(recall,            na.rm = TRUE),
    f05_score         = mean(f05_score,         na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc               = mean(mcc,               na.rm = TRUE),
    mse               = mean(mse,               na.rm = TRUE),
    rmse              = mean(rmse,              na.rm = TRUE)
  ) %>%
  ungroup()

df_avg <- df_thresh %>%
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall, f05_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold, names_to = "metric", values_to = "value")

df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

best_discrete           <- df_wide %>% slice_max(f05_score, n = 1)
best_discrete_threshold <- best_discrete$threshold
best_discrete_threshold


## ---- 5. Build result_summary ----
df_removed <- df %>%
  filter(removed == 1) %>%
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%
  mutate(predicted_bin_sigm  = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>%
  mutate(original_binary     = if_else(original_links > 0, 1, 0))

result_summary <- df_removed %>%
  group_by(train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    specificity       = TN / (TN + FP),
    precision         = TP / (TP + FP),
    recall            = TP / (TP + FN),
    f05_score         = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    balanced_accuracy = (recall + specificity) / 2,
    nse  = 1 - sum((predicted_values - original_links)^2, na.rm = TRUE) /
      sum((original_links - mean(original_links, na.rm = TRUE))^2, na.rm = TRUE),
    nnse = 1 / (2 - nse)
  ) %>%
  ungroup() %>%
  group_by(train_layer, test_layer) %>%
  summarise(
    TP                = mean(TP,                na.rm = TRUE),
    FN                = mean(FN,                na.rm = TRUE),
    TN                = mean(TN,                na.rm = TRUE),
    FP                = mean(FP,                na.rm = TRUE),
    specificity       = mean(specificity,       na.rm = TRUE),
    precision         = mean(precision,         na.rm = TRUE),
    recall            = mean(recall,            na.rm = TRUE),
    f05_score         = mean(f05_score,         na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    nse               = mean(nse,               na.rm = TRUE),
    nnse              = mean(nnse,              na.rm = TRUE)
  ) %>%
  ungroup()

result_summary <- result_summary %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Single location",
    train_layer != test_layer ~ "Added location"
  ))

result_summary <- result_summary %>%
  left_join(layer_year_df, by = c("train_layer" = "layer_num")) %>%
  rename(train_layer_name = year) %>%
  left_join(layer_year_df, by = c("test_layer" = "layer_num")) %>%
  rename(test_layer_name = year)


## ---- Fig. S20: Non-thresholded evaluation ----
# AUC-ROC and AUC-PR heatmaps -> pr_roc_host_parasite.pdf

df_eval <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1L, 0L)
  ) %>%
  group_by(train_layer, test_layer, itr) %>%
  summarise(
    auc_roc = tryCatch({
      roc_obj <- roc(response  = original_binary,
                     predictor = predicted_prob,
                     quiet = TRUE, na.rm = TRUE,
                     levels = c(0, 1), direction = "<")
      as.numeric(auc(roc_obj))
    }, error = function(e) NA_real_),
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
  group_by(train_layer, test_layer) %>%
  summarise(
    auc_roc_mean = mean(auc_roc, na.rm = TRUE),
    auc_roc_sd   = sd(auc_roc,   na.rm = TRUE),
    auc_pr_mean  = mean(auc_pr,  na.rm = TRUE),
    auc_pr_sd    = sd(auc_pr,    na.rm = TRUE),
    .groups = "drop"
  )

df_eval_summary <- df_eval_summary %>%
  left_join(layer_year_df, by = c("train_layer" = "layer_num")) %>%
  rename(train_layer_name = year) %>%
  left_join(layer_year_df, by = c("test_layer" = "layer_num")) %>%
  rename(test_layer_name = year)

lims_roc    <- range(df_eval_summary$auc_roc_mean, na.rm = TRUE)
mid_val_roc <- mean(lims_roc)

year_heatmap_auc <-
  ggplot(df_eval_summary, aes(x = train_layer_name, y = test_layer_name, fill = auc_roc_mean)) +
  geom_tile(color = "black", linewidth = 0.1) +
  geom_tile(data = df_eval_summary[df_eval_summary$train_layer == df_eval_summary$test_layer, ],
            color = "black", linewidth = 1.2) +
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "rosybrown2",
                       midpoint = mid_val_roc, na.value = "gray") +
  labs(x = "Auxiliary year", y = "Target year", fill = "ROC-AUC") +
  theme_minimal() +
  theme(
    text             = element_text(size = 18),
    plot.margin      = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed() + tme

lims_pr    <- range(df_eval_summary$auc_pr_mean, na.rm = TRUE)
mid_val_pr <- mean(lims_pr)

year_heatmap_pr <-
  ggplot(df_eval_summary, aes(x = train_layer_name, y = test_layer_name, fill = auc_pr_mean)) +
  geom_tile(color = "black", linewidth = 0.1) +
  geom_tile(data = df_eval_summary[df_eval_summary$train_layer == df_eval_summary$test_layer, ],
            color = "black", linewidth = 1.2) +
  scale_fill_gradient2(low = "lightsteelblue2", mid = "white", high = "thistle",
                       midpoint = mid_val_pr, na.value = "gray") +
  labs(x = "Auxiliary year", y = "Target year", fill = "PR-AUC") +
  theme_minimal() +
  theme(
    text             = element_text(size = 18),
    plot.margin      = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed() + tme

pr_roc <- plot_grid(
  year_heatmap_auc + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  year_heatmap_pr  + theme(axis.title.y = element_blank()),
  rel_widths = c(1, 0.92),
  labels     = c("(a)", "(b)"),
  label_size = 18,
  label_y    = 0.8
)

pdf(file   = "results/hp_analysis/paper_figs/pr_roc_host_parasite.pdf",
    width  = 13, height = 10, family = "Helvetica")
pr_roc
dev.off()


## ---- Fig. S21: Mapping potential missing links ----
# Saved as map_missing_links_host_parasite.pdf

# Overall degree per species (itr == 1 to avoid duplicates)
df_filtered <- df %>%
  filter(itr == 1, original_links != 0)

overall_host_degree <- df_filtered %>%
  group_by(node_from) %>%
  summarise(overall_host_degree = length(unique(node_to)), .groups = "drop")

overall_parasite_degree <- df_filtered %>%
  group_by(node_to) %>%
  summarise(overall_parasite_degree = length(unique(node_from)), .groups = "drop")

# Add year_pair_id and sigmoid predictions to df
df <- df %>%
  mutate(year_pair_id        = paste(train_layer, test_layer, sep = "_"),
         predicted_prob_sigm = sigmoid(predicted_values))

# Summarise observations and predictions per interaction per year-pair
df_year_pair <- df %>%
  group_by(node_from, node_to, year_pair_id) %>%
  summarise(
    observed                 = as.integer(any(original_links != 0)),
    year_pair_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )

# Summarise across all year-pairs per interaction
df_summary <- df_year_pair %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop           = mean(observed,                  na.rm = TRUE),
    avg_sigm_predicted = mean(year_pair_sigm_predicted,  na.rm = TRUE),
    n_year_pairs       = n(),
    .groups = "drop"
  )

# Proportion of YEARS (not year-pairs) in which each interaction is observed
df_filtered_self <- df_filtered %>% filter(train_layer == test_layer)

result_count <- df_filtered_self %>%
  group_by(node_to, node_from) %>%
  summarise(n_test_layers = n_distinct(test_layer), .groups = "drop") %>%
  mutate(prop_test_layers = n_test_layers / n_distinct(df$test_layer))

df_summary <- df_summary %>%
  left_join(result_count %>% select(node_to, node_from, prop_test_layers),
            by = c("node_to", "node_from")) %>%
  mutate(avg_prop_year = if_else(is.na(prop_test_layers), 0, prop_test_layers)) %>%
  select(-prop_test_layers)

# Join degree data for species ordering
df_summary <- df_summary %>%
  left_join(overall_parasite_degree, by = "node_to") %>%
  left_join(overall_host_degree,     by = "node_from")

host_order <- df_summary %>%
  distinct(node_from, overall_host_degree) %>%
  arrange(desc(overall_host_degree)) %>%
  pull(node_from)

parasite_order <- df_summary %>%
  distinct(node_to, overall_parasite_degree) %>%
  arrange(desc(overall_parasite_degree)) %>%
  pull(node_to)

df_summary$node_from <- factor(df_summary$node_from, levels = host_order)
df_summary$node_to   <- factor(df_summary$node_to,   levels = parasite_order)

map_missing_links <- ggplot(df_summary, aes(x = node_to, y = node_from)) +
  geom_tile(aes(fill = avg_prop_year)) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Observed links:\nproportion\nof years\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  new_scale_fill() +
  geom_tile(
    data  = df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
    aes(fill = avg_sigm_predicted),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0.8, 1, 0.05)) +
  theme_minimal() +
  labs(x = "Parasite", y = "Host") +
  tme +
  theme(
    axis.text.x     = element_blank(),
    axis.text.y     = element_text(size = 10),
    legend.text     = element_text(size = 12),
    legend.position = "bottom",
    legend.box      = "horizontal"
  ) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

pdf(file   = "results/hp_analysis/paper_figs/map_missing_links_host_parasite.pdf",
    width  = 11, height = 6, family = "Helvetica")
print(map_missing_links)
dev.off()


## ---- Fig. S22: Jaccard similarity and temporal distance ----
# Saved as year_jaccard_distance.pdf

# Calculate Jaccard indices for all year-pair combinations
results_jaccard <- data.frame()

for (layers_to_train in 1:num_layers) {
  for (layer_to_predict in 1:num_layers) {

    A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
    P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)

    # Jaccard parasites (rows)
    parasite_train   <- rownames(A)[rowSums(A) > 0]
    parasite_test    <- rownames(P)[rowSums(P) > 0]
    union_p          <- length(union(parasite_train, parasite_test))
    jaccard_parasite <- if (union_p == 0) NA else
      length(intersect(parasite_train, parasite_test)) / union_p

    # Jaccard hosts (columns)
    hosts_train  <- colnames(A)[colSums(A) > 0]
    hosts_test   <- colnames(P)[colSums(P) > 0]
    union_h      <- length(union(hosts_train, hosts_test))
    jaccard_hosts <- if (union_h == 0) NA else
      length(intersect(hosts_train, hosts_test)) / union_h

    # Jaccard edges
    pairs_train_strings <- apply(which(A > 0, arr.ind = TRUE), 1, function(rc)
      paste(rownames(A)[rc[1]], colnames(A)[rc[2]], sep = "_"))
    pairs_test_strings  <- apply(which(P > 0, arr.ind = TRUE), 1, function(rc)
      paste(rownames(P)[rc[1]], colnames(P)[rc[2]], sep = "_"))
    union_e       <- length(union(pairs_train_strings, pairs_test_strings))
    jaccard_edges <- if (union_e == 0) NA else
      length(intersect(pairs_train_strings, pairs_test_strings)) / union_e

    results_jaccard <- rbind(results_jaccard, data.frame(
      train_layer       = layers_to_train,
      test_layer        = layer_to_predict,
      jaccard_parasites = jaccard_parasite,
      jaccard_hosts     = jaccard_hosts,
      jaccard_edges     = jaccard_edges
    ))
  }
}

result_summary <- result_summary %>%
  left_join(results_jaccard, by = c("train_layer", "test_layer"))

# Only different-year pairs for Jaccard analysis
network_results_jaccard <- result_summary %>%
  filter(train_layer != test_layer)

# Time differences between all year pairs
time_diff_table <- expand.grid(
  train_layer = seq_along(layer_to_year),
  test_layer  = seq_along(layer_to_year),
  stringsAsFactors = FALSE
) %>%
  mutate(
    train_year      = layer_to_year[train_layer],
    test_year       = layer_to_year[test_layer],
    time_difference = abs(test_year - train_year)
  )

result_summary_temporal <- result_summary %>%
  left_join(time_diff_table %>% select(train_layer, test_layer, time_difference),
            by = c("train_layer", "test_layer"))

result_summary_temporal_dif <- result_summary_temporal %>%
  filter(train_layer != test_layer) %>%
  mutate(distance_km = time_difference)  # rename for compatibility with make_cor_plot

cor_plot_dif_year_f05 <- make_cor_plot(
  result_summary_temporal_dif, evaluator = "f05_score", extra_theme = tme
)

# Assemble Fig. S22: 2x2 grid of panels

# (a) Parasite overlap
p_a_core <- make_facet_scatter_plot2(
  data = network_results_jaccard, evaluator = "f05_score",
  pivot_cols = "jaccard_parasites", facet_scales = "free_x", tme = tme
)
label_a <- paste0("Parasite overlap: ",
                  rp_text("jaccard_parasites", "f05_score", network_results_jaccard))
p_a <- add_center_header(p_a_core, label_a, size = 12)

# (b) Host overlap (no y-axis title — right column)
p_b_core <- make_facet_scatter_plot2(
  data = network_results_jaccard, evaluator = "f05_score",
  pivot_cols = "jaccard_hosts", facet_scales = "free_x", tme = tme
) + theme(axis.title.y = element_blank())
label_b <- paste0("Host overlap: ",
                  rp_text("jaccard_hosts", "f05_score", network_results_jaccard))
p_b <- add_center_header(p_b_core, label_b, size = 12)

# (c) Edge overlap
p_c_core <- make_facet_scatter_plot2(
  data = network_results_jaccard, evaluator = "f05_score",
  pivot_cols = "jaccard_edges", facet_scales = "free_x", tme = tme
)
label_c <- paste0("Interaction overlap: ",
                  rp_text("jaccard_edges", "f05_score", network_results_jaccard))
p_c <- add_center_header(p_c_core, label_c, size = 12)

# (d) Temporal distance correlation
p_d_core <- cor_plot_dif_year_f05 +
  tme +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15))
p_d_core <- drop_text_layers(p_d_core)  # remove in-panel annotation; header carries stats
label_d  <- paste0("Temporal distance: ",
                   rp_text("distance_km", "f05_score", result_summary_temporal_dif))
p_d <- add_center_header(p_d_core, label_d, size = 12)

year_jaccard_distance <- plot_grid(
  p_a, p_b, p_c, p_d,
  labels     = c("(a)", "(b)", "(c)", "(d)"),
  ncol       = 2,
  align      = "hv",
  axis       = "tblr",
  label_size = 13,
  label_y    = c(0.96, 0.96, 0.96, 0.96),
  label_x    = c(0, 0, 0, 0),
  rel_widths = c(1, 0.95, 1, 0.95)
)

pdf(file   = "results/hp_analysis/paper_figs/year_jaccard_distance.pdf",
    width  = 8, height = 8, family = "Helvetica")
print(year_jaccard_distance)
dev.off()
