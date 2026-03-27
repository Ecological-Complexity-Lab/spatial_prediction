library(dplyr)
library(jsonlite)

# ── Load data ──────────────────────────────────────────────────────────────────
df_combined <- readRDS("results/df_combined.rds")

cat("Columns:", paste(names(df_combined), collapse=", "), "\n")
cat("node_from is factor:", is.factor(df_combined$node_from), "\n")
cat("node_to   is factor:", is.factor(df_combined$node_to),   "\n")

# ── Extract correct ordering ───────────────────────────────────────────────────
# The R analysis (lines 1663-1671) orders by desc(overall_degree).
# overall_plant_degree = number of unique pollinators per plant (from df_filtered)
# overall_poll_degree  = number of unique plants per pollinator (from df_filtered)
# df_filtered = df %>% filter(itr == 1, original_links != 0)
#
# df_combined was created AFTER plant_order/poll_order were computed and factor
# levels were set on df_summary (lines 1703-1704). So if node_from/node_to are
# factors in df_combined.rds, their levels ARE the correct order.

if (is.factor(df_combined$node_from) && is.factor(df_combined$node_to)) {
  plant_order_correct <- levels(df_combined$node_from)
  poll_order_correct  <- levels(df_combined$node_to)
  cat("\nPlant order from factor levels (n =", length(plant_order_correct), "):\n")
  print(plant_order_correct)
  cat("\nPoll order from factor levels (first 20, n =", length(poll_order_correct), "):\n")
  print(head(poll_order_correct, 20))
} else {
  cat("\nNot factors — recomputing from scratch...\n")
  # Replicate lines 1456-1485 and 1663-1671 of the analysis script
  # We need the original df for this. df_combined only has aggregated data.
  stop("node_from/node_to are not factors in df_combined.rds. Need original df.")
}

# ── Compare with current JSON ──────────────────────────────────────────────────
json_data <- fromJSON("docs/figure3a_data.json")

cat("\n--- Comparing plant_order ---\n")
json_plants <- json_data$plant_order
cat("JSON plants (n =", length(json_plants), "):\n")
print(json_plants)
cat("Correct plants (n =", length(plant_order_correct), "):\n")
print(plant_order_correct)

if (identical(json_plants, plant_order_correct)) {
  cat("plant_order: MATCH\n")
} else {
  cat("plant_order: MISMATCH\n")
  diffs <- which(json_plants != plant_order_correct)
  cat("First mismatches at positions:", head(diffs, 10), "\n")
  for (i in head(diffs, 5)) {
    cat(sprintf("  pos %d: JSON='%s'  correct='%s'\n", i, json_plants[i], plant_order_correct[i]))
  }
}

cat("\n--- Comparing poll_order (first 20) ---\n")
json_polls <- json_data$poll_order
cat("JSON polls (n =", length(json_polls), "):\n")
print(head(json_polls, 20))
cat("Correct polls (n =", length(poll_order_correct), "):\n")
print(head(poll_order_correct, 20))

if (identical(json_polls, poll_order_correct)) {
  cat("poll_order: MATCH\n")
} else {
  cat("poll_order: MISMATCH\n")
  diffs <- which(json_polls != poll_order_correct)
  cat("First mismatches at positions:", head(diffs, 10), "\n")
  for (i in head(diffs, 10)) {
    cat(sprintf("  pos %d: JSON='%s'  correct='%s'\n", i, json_polls[i], poll_order_correct[i]))
  }
}
