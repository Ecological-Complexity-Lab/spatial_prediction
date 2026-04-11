# Export data for interactive Figure 4
# Generates: interactive_figure4/figure4_data.json
#
# Pipeline copied verbatim from Abramov_et_al_spatial_prediction_analysis.R,
# with thresholds expanded to a finer grid for the interactive slider.
# Only change vs. the paper: thresholds <- seq(0.50, 0.90, by = 0.025)
#                            instead of  seq(0, 1, by = 0.1)

library(dplyr)
library(tidyr)
library(jsonlite)
library(PRROC)

sigmoid <- function(x) 1 / (1 + exp(-x))

layer_names <- c(
  "1" = "Tenerife Teno",  "2" = "Tenerife South",
  "3" = "Gran Canaria",   "4" = "Fuerteventura",
  "5" = "Western Sahara", "6" = "El Hierro",
  "7" = "La Gomera"
)

# ---- Load & filter (lines 726-736 of paper script) ----
cat("Loading raw predictions...\n")
combined_results <- readRDS("results/predictions_island_scale.rds")
cat("Loaded", nrow(combined_results), "rows\n")

combined_results <- combined_results %>%
  filter(k == 2) %>%
  filter(!(input_lambda %in% c(1, 5, 50, 100)))

df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

rm(combined_results); gc()

# ---- Threshold analysis (lines 748-798 of paper script) ----
# Only change: finer threshold grid for the interactive slider
thresholds <- seq(0.50, 0.90, by = 0.025)

df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1, 0)
  )

cat("Computing metrics across", length(thresholds), "thresholds...\n")
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP          = sum(original_binary == 1 & predicted_bin == 1),
    FN          = sum(original_binary == 1 & predicted_bin == 0),
    TN          = sum(original_binary == 0 & predicted_bin == 0),
    FP          = sum(original_binary == 0 & predicted_bin == 1),
    specificity = TN / (TN + FP),
    precision   = TP / (TP + FP),
    recall      = TP / (TP + FN),
    f05_score   = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    f1_score    = (2    ) * (precision * recall) / (precision + recall),
    .groups = "drop"
  ) %>%
  group_by(emln_id, train_layer, test_layer, threshold) %>%
  summarise(
    f05_score = mean(f05_score, na.rm = TRUE),
    f1_score  = mean(f1_score,  na.rm = TRUE),
    .groups = "drop"
  )

cat("df_thresh rows:", nrow(df_thresh), "\n")

# ---- PR-AUC (lines 841-877 of paper script) ----
cat("Computing PR-AUC...\n")
df_eval <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob  = sigmoid(predicted_values),
    original_binary = if_else(original_links > 0, 1L, 0L)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    auc_pr = tryCatch({
      pos <- predicted_prob[original_binary == 1]
      neg <- predicted_prob[original_binary == 0]
      if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
      pr_obj <- pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)
      pr_obj$auc.integral
    }, error = function(e) NA_real_),
    .groups = "drop"
  )

df_eval_summary <- df_eval %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    auc_pr_mean = mean(auc_pr, na.rm = TRUE),
    .groups = "drop"
  )

rm(df_prepped, df_eval, df); gc()

# ---- Join metadata (Jaccard, distance, layer names) ----
meta <- read.csv("results/result_summary_island.csv") %>%
  select(train_layer, test_layer,
         train_layer_name, test_layer_name,
         jaccard_pollinators, jaccard_plants, jaccard_edges,
         distance_km, layer_comparison) %>%
  distinct()

df_full <- df_thresh %>%
  left_join(meta,           by = c("train_layer", "test_layer")) %>%
  left_join(df_eval_summary, by = c("emln_id", "train_layer", "test_layer")) %>%
  arrange(emln_id, train_layer, test_layer, threshold)

# ---- Build JSON ----
# scatter: one entry per (emln_id, train_layer, test_layer) — matches paper's data points
# pairs:   averaged across emln_id — for the heatmap
cat("Building JSON...\n")

scatter_keys <- df_full %>%
  select(emln_id, train_layer, test_layer) %>%
  distinct() %>%
  arrange(emln_id, train_layer, test_layer)

scatter_list <- vector("list", nrow(scatter_keys))
for (i in seq_len(nrow(scatter_keys))) {
  eid <- scatter_keys$emln_id[i]
  tl  <- scatter_keys$train_layer[i]
  tst <- scatter_keys$test_layer[i]

  rows     <- df_full %>% filter(emln_id == eid, train_layer == tl, test_layer == tst) %>% arrange(threshold)
  meta_row <- rows[1, ]

  metrics_by_thresh <- setNames(
    lapply(seq_len(nrow(rows)), function(j)
      list(f05 = round(rows$f05_score[j], 6), f1 = round(rows$f1_score[j], 6))
    ),
    as.character(round(rows$threshold, 3))
  )

  scatter_list[[i]] <- list(
    emln_id     = eid,
    train_layer = tl,
    test_layer  = tst,
    train_name  = layer_names[[as.character(tl)]],
    test_name   = layer_names[[as.character(tst)]],
    layer_comparison    = meta_row$layer_comparison,
    jaccard_pollinators = if (!is.na(meta_row$jaccard_pollinators)) round(meta_row$jaccard_pollinators, 6) else NULL,
    jaccard_plants      = if (!is.na(meta_row$jaccard_plants))      round(meta_row$jaccard_plants, 6)      else NULL,
    jaccard_edges       = if (!is.na(meta_row$jaccard_edges))       round(meta_row$jaccard_edges, 6)       else NULL,
    distance_km         = round(meta_row$distance_km, 3),
    pr_auc              = if (!is.na(meta_row$auc_pr_mean))         round(meta_row$auc_pr_mean, 6)         else NULL,
    metrics = metrics_by_thresh
  )
}

# Heatmap: average across emln_id (same as scatter if one lambda per pair)
df_pair <- df_full %>%
  group_by(train_layer, test_layer, threshold) %>%
  summarise(f05_score = mean(f05_score, na.rm = TRUE),
            f1_score  = mean(f1_score,  na.rm = TRUE), .groups = "drop")

df_pair_meta <- df_full %>%
  group_by(train_layer, test_layer) %>%
  summarise(auc_pr_mean = mean(auc_pr_mean, na.rm = TRUE), .groups = "drop") %>%
  left_join(meta, by = c("train_layer", "test_layer"))

all_pairs <- df_pair %>% select(train_layer, test_layer) %>% distinct() %>% arrange(train_layer, test_layer)
pairs_list <- vector("list", nrow(all_pairs))
for (i in seq_len(nrow(all_pairs))) {
  tl  <- all_pairs$train_layer[i]
  tst <- all_pairs$test_layer[i]

  rows     <- df_pair %>% filter(train_layer == tl, test_layer == tst) %>% arrange(threshold)
  meta_row <- df_pair_meta %>% filter(train_layer == tl, test_layer == tst)

  metrics_by_thresh <- setNames(
    lapply(seq_len(nrow(rows)), function(j)
      list(f05 = round(rows$f05_score[j], 6), f1 = round(rows$f1_score[j], 6))
    ),
    as.character(round(rows$threshold, 3))
  )

  pairs_list[[i]] <- list(
    train_layer = tl,
    test_layer  = tst,
    train_name  = layer_names[[as.character(tl)]],
    test_name   = layer_names[[as.character(tst)]],
    layer_comparison    = meta_row$layer_comparison,
    jaccard_pollinators = if (!is.na(meta_row$jaccard_pollinators)) round(meta_row$jaccard_pollinators, 6) else NULL,
    jaccard_plants      = if (!is.na(meta_row$jaccard_plants))      round(meta_row$jaccard_plants, 6)      else NULL,
    jaccard_edges       = if (!is.na(meta_row$jaccard_edges))       round(meta_row$jaccard_edges, 6)       else NULL,
    distance_km         = round(meta_row$distance_km, 3),
    pr_auc              = if (!is.na(meta_row$auc_pr_mean))         round(meta_row$auc_pr_mean, 6)         else NULL,
    metrics = metrics_by_thresh
  )
}

output <- list(
  thresholds        = round(thresholds, 3),
  layer_names       = as.list(layer_names),
  default_threshold = 0.7,
  pairs             = pairs_list,
  scatter           = scatter_list
)

json_path <- "interactive_figure4/figure4_data.json"
write(toJSON(output, auto_unbox = TRUE, pretty = FALSE), json_path)
cat("Written to", json_path, "\n")
cat("File size:", round(file.size(json_path) / 1024, 1), "KB\n")
cat("Scatter entries:", length(scatter_list), " Pairs:", length(pairs_list), "\n")
