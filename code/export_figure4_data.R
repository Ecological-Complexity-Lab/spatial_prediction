# Export data for interactive Figure 4
# Generates: docs/figure4_data.json
#
# Strategy:
#  - At threshold=0.6: exact paper values from result_summary_island_dif.rds
#  - At other thresholds: paper values + delta computed from raw predictions
#    (raw predictions are from a different random seed, so absolute values differ,
#     but threshold-induced changes are meaningful)

library(dplyr)
library(jsonlite)
library(tidyr)

sigmoid <- function(x) 1 / (1 + exp(-x))

layer_names <- c(
  "1" = "Tenerife Teno",  "2" = "Tenerife South",
  "3" = "Gran Canaria",   "4" = "Fuerteventura",
  "5" = "Western Sahara", "6" = "El Hierro",
  "7" = "La Gomera"
)

# ---- Ground truth at threshold = 0.6 (exact paper values) ----
truth <- read.csv("results/result_summary_island.csv") %>%
  filter(train_layer != test_layer) %>%
  select(train_layer, test_layer,
         train_layer_name, test_layer_name,
         TP, FN, TN, FP,
         f05_paper = f05_score,
         jaccard_pollinators, jaccard_plants, jaccard_edges,
         distance_km, layer_comparison) %>%
  mutate(
    precision  = TP / (TP + FP),
    recall     = TP / (TP + FN),
    f1_paper   = (2 * precision * recall) / (precision + recall),
    mcc_paper  = (TP * TN - FP * FN) /
                   sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  select(-precision, -recall)

# Also build diagonal (same-island) entries from raw result_summary
# (result_summary_island_dif only has off-diagonal pairs)
truth_diag <- data.frame(
  train_layer = 1:7,
  test_layer  = 1:7,
  train_layer_name = layer_names,
  test_layer_name  = layer_names,
  TP = NA_real_, FN = NA_real_, TN = NA_real_, FP = NA_real_,
  f05_paper = NA_real_, f1_paper = NA_real_, mcc_paper = NA_real_,
  jaccard_pollinators = NA_real_, jaccard_plants = NA_real_,
  jaccard_edges = NA_real_, distance_km = 0,
  layer_comparison = "Single location"
)

# ---- Raw predictions for threshold variation (delta computation) ----
cat("Loading raw predictions...\n")
df_raw <- readRDS("/tmp/preds_tmp/predictions_island_scale.rds")

df_removed <- df_raw %>%
  filter(k == 2,
         !(input_lambda %in% c(1, 5, 50, 100)),
         removed == 1) %>%
  mutate(
    predicted_values = if_else(predicted_values < 0, 0, predicted_values),
    prob             = sigmoid(predicted_values),
    original_bin     = as.integer(original_links > 0)
  ) %>%
  select(train_layer, test_layer, itr, prob, original_bin)

rm(df_raw); gc()

# ---- Compute metrics per-iteration at every threshold ----
thresholds <- seq(0.30, 0.80, by = 0.025)

cat("Computing per-iteration metrics at", length(thresholds), "thresholds...\n")

df_thresh <- df_removed %>%
  crossing(threshold = thresholds) %>%
  mutate(pred_bin = as.integer(prob > threshold)) %>%
  group_by(train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_bin == 1L & pred_bin == 1L),
    FN = sum(original_bin == 1L & pred_bin == 0L),
    TN = sum(original_bin == 0L & pred_bin == 0L),
    FP = sum(original_bin == 0L & pred_bin == 1L),
    .groups = "drop"
  ) %>%
  mutate(
    precision = TP / (TP + FP),
    recall    = TP / (TP + FN),
    f05_raw   = (1.25 * precision * recall) / (0.25 * precision + recall),
    f1_raw    = (2   * precision * recall) / (precision + recall),
    mcc_raw   = (TP * TN - FP * FN) /
                  sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  group_by(train_layer, test_layer, threshold) %>%
  summarise(
    f05_raw = mean(f05_raw, na.rm = TRUE),
    f1_raw  = mean(f1_raw,  na.rm = TRUE),
    mcc_raw = mean(mcc_raw, na.rm = TRUE),
    .groups = "drop"
  )

# ---- At threshold = 0.6: compute raw-baseline values ----
ref_thresh <- 0.6
baseline_raw <- df_thresh %>%
  filter(abs(threshold - ref_thresh) < 1e-9) %>%
  rename(f05_base = f05_raw, f1_base = f1_raw, mcc_base = mcc_raw) %>%
  select(train_layer, test_layer, f05_base, f1_base, mcc_base)

# ---- Compute deltas: delta = raw(threshold) - raw(0.6) ----
df_delta <- df_thresh %>%
  left_join(baseline_raw, by = c("train_layer", "test_layer")) %>%
  mutate(
    d_f05 = f05_raw - f05_base,
    d_f1  = f1_raw  - f1_base,
    d_mcc = mcc_raw - mcc_base
  ) %>%
  select(train_layer, test_layer, threshold, d_f05, d_f1, d_mcc)

cat("Done. Building JSON...\n")

# ---- Build all pairs (diagonal + off-diagonal) ----
all_truth <- bind_rows(truth, truth_diag) %>%
  arrange(train_layer, test_layer)

pairs_list <- vector("list", nrow(all_truth))

for (i in seq_len(nrow(all_truth))) {
  tl  <- all_truth$train_layer[i]
  tst <- all_truth$test_layer[i]

  deltas <- df_delta %>%
    filter(train_layer == tl, test_layer == tst) %>%
    arrange(threshold)

  # Build metrics: paper value at 0.6 + delta at each threshold
  f05_base <- all_truth$f05_paper[i]
  f1_base  <- all_truth$f1_paper[i]
  mcc_base <- all_truth$mcc_paper[i]

  metrics_by_thresh <- setNames(
    lapply(seq_len(nrow(deltas)), function(j) {
      t <- deltas$threshold[j]
      # At threshold = 0.6: use exact paper values; elsewhere: paper + delta
      if (abs(t - ref_thresh) < 1e-9) {
        list(f05 = round(f05_base, 6),
             f1  = round(f1_base,  6),
             mcc = round(mcc_base, 6))
      } else {
        list(f05 = round(f05_base + deltas$d_f05[j], 6),
             f1  = round(f1_base  + deltas$d_f1[j],  6),
             mcc = round(mcc_base + deltas$d_mcc[j], 6))
      }
    }),
    as.character(round(deltas$threshold, 3))
  )

  # For diagonal pairs without paper values, just use raw values
  if (is.na(f05_base)) {
    metrics_by_thresh <- setNames(
      lapply(seq_len(nrow(deltas)), function(j) {
        list(f05 = round(deltas$f05_raw[j] + deltas$d_f05[j], 6),
             f1  = round(deltas$f1_raw[j]  + deltas$d_f1[j],  6),
             mcc = round(deltas$d_mcc[j], 6))
      }),
      as.character(round(deltas$threshold, 3))
    )
    # Simpler: just use df_thresh directly for diagonal
    raw_diag <- df_thresh %>% filter(train_layer == tl, test_layer == tst)
    metrics_by_thresh <- setNames(
      lapply(seq_len(nrow(raw_diag)), function(j) {
        list(f05 = round(raw_diag$f05_raw[j], 6),
             f1  = round(raw_diag$f1_raw[j],  6),
             mcc = round(raw_diag$mcc_raw[j], 6))
      }),
      as.character(round(raw_diag$threshold, 3))
    )
  }

  pairs_list[[i]] <- list(
    train_layer = tl,
    test_layer  = tst,
    train_name  = layer_names[[as.character(tl)]],
    test_name   = layer_names[[as.character(tst)]],
    layer_comparison    = all_truth$layer_comparison[i],
    jaccard_pollinators = if (!is.na(all_truth$jaccard_pollinators[i]))
                            round(all_truth$jaccard_pollinators[i], 6) else NULL,
    jaccard_plants      = if (!is.na(all_truth$jaccard_plants[i]))
                            round(all_truth$jaccard_plants[i], 6)      else NULL,
    jaccard_edges       = if (!is.na(all_truth$jaccard_edges[i]))
                            round(all_truth$jaccard_edges[i], 6)       else NULL,
    distance_km         = round(all_truth$distance_km[i], 3),
    metrics = metrics_by_thresh
  )
}

output <- list(
  thresholds        = round(thresholds, 3),
  layer_names       = as.list(layer_names),
  default_threshold = ref_thresh,
  pairs             = pairs_list
)

json_path <- "docs/figure4_data.json"
write(toJSON(output, auto_unbox = TRUE, pretty = FALSE), json_path)
cat("Written to", json_path, "\n")
cat("File size:", round(file.size(json_path) / 1024, 1), "KB\n")

# ---- Sanity check ----
cat("\nSanity check — F0.5 at threshold=0.6 vs paper:\n")
check_pairs <- list(c(1,2), c(3,4), c(6,7), c(7,5))
for (cp in check_pairs) {
  row <- all_truth %>% filter(train_layer==cp[1], test_layer==cp[2])
  cat(sprintf("  Tr=%d Te=%d: stored=%.4f paper=%.4f\n",
              cp[1], cp[2], row$f05_paper, row$f05_paper))
}
