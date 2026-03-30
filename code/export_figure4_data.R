# Export data for interactive Figure 4
# Generates: results/result_summary_all_thresholds.csv
#            docs/figure4_data.json
#
# Pipeline: exact replication of paper analysis
#   filter k==2, exclude extreme lambdas, removed==1,
#   clamp negatives, sigmoid, per-iteration confusion matrix,
#   compute metrics per iteration, average across iterations.
# Jaccard + distance metadata joined from result_summary_island.csv.

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

# ---- Load raw predictions ----
cat("Loading raw predictions...\n")
zip_path <- "results/predictions_island_scale.zip"
tmp_dir  <- tempdir()
unzip(zip_path, exdir = tmp_dir)
rds_file <- list.files(tmp_dir, pattern = "\\.rds$", full.names = TRUE)[1]
df_raw   <- readRDS(rds_file)
cat("Loaded", nrow(df_raw), "rows\n")

# ---- Apply paper pipeline ----
df_removed <- df_raw %>%
  filter(k == 2,
         !(input_lambda %in% c(1, 5, 50, 100)),
         removed == 1) %>%
  mutate(
    predicted_values = pmax(predicted_values, 0),
    prob             = sigmoid(predicted_values),
    original_bin     = as.integer(original_links > 0)
  ) %>%
  select(train_layer, test_layer, itr, prob, original_bin)

rm(df_raw); gc()
cat("After filtering:", nrow(df_removed), "rows\n")

# ---- Compute metrics at every threshold ----
thresholds <- seq(0.30, 0.80, by = 0.025)
cat("Computing metrics at", length(thresholds), "thresholds...\n")

df_all <- df_removed %>%
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
    f05  = (1.25 * precision * recall) / (0.25 * precision + recall),
    f1   = (2   * precision * recall) / (precision + recall),
    mcc  = (TP * TN - FP * FN) /
             sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  group_by(train_layer, test_layer, threshold) %>%
  summarise(
    f05 = mean(f05, na.rm = TRUE),
    f1  = mean(f1,  na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    .groups = "drop"
  )

cat("Metrics computed:", nrow(df_all), "rows\n")

# ---- Join Jaccard + distance metadata ----
meta <- read.csv("results/result_summary_island.csv") %>%
  select(train_layer, test_layer,
         train_layer_name, test_layer_name,
         jaccard_pollinators, jaccard_plants, jaccard_edges,
         distance_km, layer_comparison) %>%
  distinct()

df_final <- df_all %>%
  left_join(meta, by = c("train_layer", "test_layer")) %>%
  arrange(train_layer, test_layer, threshold)

# ---- Save CSV ----
csv_path <- "results/result_summary_all_thresholds.csv"
write.csv(df_final, csv_path, row.names = FALSE)
cat("Written to", csv_path, "\n")
cat("Rows:", nrow(df_final), "  Pairs:", n_distinct(df_final[c("train_layer","test_layer")]), "\n")

# ---- Build JSON ----
cat("Building JSON...\n")

all_pairs <- df_final %>%
  select(train_layer, test_layer) %>%
  distinct() %>%
  arrange(train_layer, test_layer)

pairs_list <- vector("list", nrow(all_pairs))

for (i in seq_len(nrow(all_pairs))) {
  tl  <- all_pairs$train_layer[i]
  tst <- all_pairs$test_layer[i]

  rows <- df_final %>%
    filter(train_layer == tl, test_layer == tst) %>%
    arrange(threshold)

  meta_row <- rows[1, ]

  metrics_by_thresh <- setNames(
    lapply(seq_len(nrow(rows)), function(j) {
      list(f05 = round(rows$f05[j], 6),
           f1  = round(rows$f1[j],  6),
           mcc = round(rows$mcc[j], 6))
    }),
    as.character(round(rows$threshold, 3))
  )

  pairs_list[[i]] <- list(
    train_layer = tl,
    test_layer  = tst,
    train_name  = layer_names[[as.character(tl)]],
    test_name   = layer_names[[as.character(tst)]],
    layer_comparison    = meta_row$layer_comparison,
    jaccard_pollinators = if (!is.na(meta_row$jaccard_pollinators))
                            round(meta_row$jaccard_pollinators, 6) else NULL,
    jaccard_plants      = if (!is.na(meta_row$jaccard_plants))
                            round(meta_row$jaccard_plants, 6)      else NULL,
    jaccard_edges       = if (!is.na(meta_row$jaccard_edges))
                            round(meta_row$jaccard_edges, 6)       else NULL,
    distance_km         = round(meta_row$distance_km, 3),
    metrics = metrics_by_thresh
  )
}

output <- list(
  thresholds        = round(thresholds, 3),
  layer_names       = as.list(layer_names),
  default_threshold = 0.7,
  pairs             = pairs_list
)

json_path <- "docs/figure4_data.json"
write(toJSON(output, auto_unbox = TRUE, pretty = FALSE), json_path)
cat("Written to", json_path, "\n")
cat("File size:", round(file.size(json_path) / 1024, 1), "KB\n")

# ---- Sanity check ----
cat("\nSanity check — F0.5 at threshold=0.7 for selected pairs:\n")
check_pairs <- list(c(1,2), c(3,4), c(6,7), c(7,5))
for (cp in check_pairs) {
  row <- df_final %>% filter(train_layer==cp[1], test_layer==cp[2],
                              abs(threshold - 0.7) < 1e-9)
  cat(sprintf("  Tr=%d Te=%d: f05=%.4f  f1=%.4f  mcc=%.4f\n",
              cp[1], cp[2], row$f05, row$f1, row$mcc))
}
