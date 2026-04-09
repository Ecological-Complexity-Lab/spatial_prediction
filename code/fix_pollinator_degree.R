library(dplyr)
library(jsonlite)

# ── Load and filter predictions — replicating main script exactly ──────────
# Main script: results_file <- "results/predictions_island_scale.rds" (line 580)
# Saved BEFORE k/lambda filtering, so we apply the same filters here.
# combined_results filters (lines 727-729):
#   filter(k == 2)
#   filter(!(input_lambda %in% c(1, 5, 50, 100)))
# df_filtered (line 1447-1448):
#   filter(itr == 1, original_links != 0)

df_raw <- readRDS("results/predictions_island_scale.rds")

overall_poll_degree <- df_raw %>%
  filter(k == 2,
         !(input_lambda %in% c(1, 5, 50, 100)),
         itr == 1,
         original_links != 0) %>%
  group_by(node_to) %>%
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

cat("Pollinators with degree computed:", nrow(overall_poll_degree), "\n")
rm(df_raw); gc()

# ── Patch global_degree for pollinators in figure3a_data.json ─────────────
json_fig3 <- fromJSON("docs/figure3a_data.json", simplifyVector = FALSE)

n_updated <- 0
for (i in seq_along(json_fig3$pollinators)) {
  sp <- json_fig3$pollinators[[i]]$id
  correct_deg <- overall_poll_degree$overall_poll_degree[overall_poll_degree$node_to == sp]
  if (length(correct_deg) == 1) {
    json_fig3$pollinators[[i]]$global_degree <- correct_deg
    n_updated <- n_updated + 1
  }
}

write(toJSON(json_fig3, auto_unbox = TRUE, null = "null", digits = 4),
      "docs/figure3a_data.json")

cat("Updated global_degree for", n_updated, "of", length(json_fig3$pollinators),
    "pollinators in docs/figure3a_data.json\n")
