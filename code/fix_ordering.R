library(dplyr)
library(jsonlite)

# ── Load data ──────────────────────────────────────────────────────────────────
df_combined <- readRDS("results/df_combined.rds")

# The factor levels in df_combined ARE the correct plot order:
# df_summary$node_from and node_to had their levels set from plant_order/poll_order
# (lines 1703-1704 of the analysis), which were computed by arrange(desc(degree)).
plant_order_correct <- levels(df_combined$node_from)
poll_order_correct  <- levels(df_combined$node_to)

# ── Load current JSON ──────────────────────────────────────────────────────────
json_data <- fromJSON("docs/figure3a_data.json", simplifyVector = FALSE)

# ── Replace orderings ──────────────────────────────────────────────────────────
json_data$plant_order <- as.list(plant_order_correct)
json_data$poll_order  <- as.list(poll_order_correct)

# Also reorder the plants/pollinators arrays to match (cosmetic, but consistent)
plant_id_order <- match(plant_order_correct,
                        sapply(json_data$plants, function(x) x$id))
json_data$plants <- json_data$plants[plant_id_order]

poll_id_order <- match(poll_order_correct,
                       sapply(json_data$pollinators, function(x) x$id))
json_data$pollinators <- json_data$pollinators[poll_id_order]

# ── Write updated JSON ─────────────────────────────────────────────────────────
write(toJSON(json_data, auto_unbox = TRUE, null = "null", digits = 4),
      "docs/figure3a_data.json")

cat("Done. plant_order and poll_order updated.\n")
cat("Plants:", length(plant_order_correct), "\n")
cat("Pollinators:", length(poll_order_correct), "\n")

# Quick sanity check
json_check <- fromJSON("docs/figure3a_data.json")
cat("JSON plant_order[1]:", json_check$plant_order[1], "\n")  # should be Euphorbia_balsamifera_m
cat("JSON plant_order[39]:", json_check$plant_order[39], "\n") # should be Neochamaelea_pulverulenta
cat("JSON poll_order[15]:", json_check$poll_order[15], "\n")   # should be Chrysomya_albiceps
