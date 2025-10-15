library(dplyr)
library(tidyr)

# ---- configurable threshold + helper ----
sigmoid <- function(x) 1/(1 + exp(-x))
best_optimal_threshold <- 0.6
# =========================================================
# 1) WITHIN-LAYER: all unique existing links (original_links != 0)
# =========================================================

# check if nodes are the same in self prediction and added location
layer1 <- df %>% filter(test_layer == 1 & itr == 1 & train_layer == 2)
view(layer1)
layer1_self <- df %>% filter(test_layer == 1 & itr == 1 & train_layer == 1)
view(layer1_self)

layer1_polls <- unique(layer1$node_to)
layer1_self_polls <- unique(layer1_self$node_to)

intersect(layer1_polls, layer1_self_polls)
length(intersect(layer1_polls, layer1_self_polls))
length(layer1_polls) # they're the same. so we can just compare link prediction between cross-layer combos

layer1_plants <- unique(layer1$node_from)
layer1_self_plants <- unique(layer1_self$node_from)

intersect(layer1_plants, layer1_self_plants)
length(intersect(layer1_plants, layer1_self_plants))
length(layer1_plants) # extra check for plants


# 
# within_existing <- df %>%
#   filter(train_layer == test_layer) %>%
#   transmute(layer = train_layer, node_from, node_to,
#             present = original_links != 0) %>%
#   group_by(layer, node_from, node_to) %>%
#   summarise(present = any(present, na.rm = TRUE), .groups = "drop") %>%
#   filter(present)   # keep only existing links

# =========================================================
# 2) CROSS-LAYER: tag test-layer EXISTING links as unique/shared
#    (unique = present in test only; shared = present in test AND train)
# =========================================================
# present_test  <- within_existing %>% rename(test_layer  = layer, in_test  = present)
# present_train <- within_existing %>% rename(train_layer = layer, in_train = present)

library(dplyr)
library(tidyr)

# 1) For each layer (using test_layer as the layer ID), mark which pairs exist (original_links != 0)
links_by_layer <- df %>%
  group_by(test_layer, node_from, node_to) %>%
  summarise(is_link = any(original_links != 0, na.rm = TRUE), .groups = "drop")

# 2) Join those "existence" flags for the current TEST layer (P) and the TRAIN layer (A)
df_labeled <- df %>%
  # existence of this pair in the TEST layer P
  left_join(links_by_layer %>%
              rename(in_P = is_link),
            by = c("test_layer", "node_from", "node_to")) %>%
  # existence of this pair in the TRAIN layer A (note: compare via the other test layer's set)
  left_join(links_by_layer %>%
              rename(train_layer = test_layer, in_A = is_link),
            by = c("train_layer", "node_from", "node_to")) %>%
  mutate(
    in_P = coalesce(in_P, FALSE),
    in_A = coalesce(in_A, FALSE),
    overlap_label = case_when(
      original_links == 0               ~ "non-link",     # absent in current TEST layer
      in_P & in_A                       ~ "shared",
      in_P & !in_A                      ~ "unique_to_P",
      !in_P & in_A                      ~ "unique_to_A",
      TRUE                              ~ "non-link"      # safety fallback
    )
  )

# focus only on cross-layer comparisons
df_labeled_cross <- df_labeled %>% filter(train_layer != test_layer)

# (Optional) Quick sanity counts
label_counts <- df_labeled_cross %>%
  count(train_layer, test_layer, overlap_label, name = "n")

label_counts

# evaluate
cross_eval <- df_labeled_cross %>%
  filter(removed == 1) %>%
  mutate(
    y_true = as.integer(original_links != 0),
    y_prob = sigmoid(predicted_values),
    y_pred = as.integer(y_prob >= best_optimal_threshold)
  )


# Shows counts of 0/1 by category; negatives will appear under NA (expected)
with(cross_eval, table(overlap_label, y_true, useNA = "ifany")) # we have about 10 times more unique links than shared

metric_fun <- function(d) {
  tp <- sum(d$y_pred == 1 & d$y_true == 1, na.rm = TRUE)
  fp <- sum(d$y_pred == 1 & d$y_true == 0, na.rm = TRUE)
  fn <- sum(d$y_pred == 0 & d$y_true == 1, na.rm = TRUE)
  recall <- (tp / (tp + fn))
  precision <- tp / (tp + fp)
  f1 <- 2 * (precision * recall) / (precision + recall)
  tibble(tp = tp, fp = fp, fn = fn, f1 = f1)
}

# F1 for SHARED: use all removed negatives + removed positives tagged "shared"
res_shared <- cross_eval %>%
  filter(y_true == 0 | overlap_label == "shared") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(overlap_label = "shared")

# F1 for UNIQUE: use all removed negatives + removed positives tagged "unique"
res_unique <- cross_eval %>%
  filter(y_true == 0 | overlap_label == "unique_to_P") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(overlap_label = "unique_to_P")

f1_by_itr <- bind_rows(res_shared, res_unique)

f1_by_itr %>% 
  mutate(predicted_positives = tp + fp) %>%    # l
  group_by(overlap_label) %>% 
  summarise(
    mean_true_positives = mean(tp),            # m
    count_observations = n(),                  # n
    total_predicted_positives = sum(predicted_positives)  # s
  )
# ---- compare ----
df_plot_f1 <- f1_by_itr %>%
  filter(overlap_label %in% c("unique_to_P","shared"), !is.na(f1)) %>%
  mutate(overlap_label = factor(overlap_label, levels = c("unique_to_P","shared")))

# --- Basic boxplot across all iters/combos ---
shared_vs_unique_links <- ggplot(df_plot_f1, aes(x = overlap_label, y = f1, fill = overlap_label)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    notch = TRUE,
    color = "grey30" # border color of boxplots
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.2,
    size = 1.6,
    color = "lightsteelblue"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 2.5,
    shape = 23,
    fill = "white"
  ) +
  labs(x = NULL, y = "F1 score") +
  stat_compare_means(
    method        = "wilcox.test",
    label         = "p.format",
    p.format.args = list(digits = 2, scientific = TRUE),
    label.y       = Inf,
    vjust         = 1.5,
    label.x       = 1.5,
    tip.length    = 0.01,
    size          = 3.5
  ) +
  scale_fill_manual(
    values = c("unique_to_P" = "lightsteelblue2",
               "shared"      = "wheat2")
  ) +
  scale_x_discrete(
    labels = c("Unique to P", "Shared")
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + tme

# pdf(
#   file   = "shared_vs_unique_links.pdf",
#   width  = 4,
#   height = 4,
#   family = "Helvetica"
# )
# print(shared_vs_unique_links)
# dev.off()     # close the file
