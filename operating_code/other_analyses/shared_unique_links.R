best_optimal_threshold <- 0.6
library(dplyr)
library(tidyr)

# ---- configurable threshold + helper ----
sigmoid <- function(x) 1/(1 + exp(-x))

# =========================================================
# 1) WITHIN-LAYER: all unique existing links (original_links != 0)
# =========================================================
within_existing <- df %>%
  filter(train_layer == test_layer) %>%
  transmute(layer = train_layer, node_from, node_to,
            present = original_links != 0) %>%
  group_by(layer, node_from, node_to) %>%
  summarise(present = any(present, na.rm = TRUE), .groups = "drop") %>%
  filter(present)   # keep only existing links

# =========================================================
# 2) CROSS-LAYER: tag test-layer EXISTING links as unique/shared
#    (unique = present in test only; shared = present in test AND train)
# =========================================================
present_test  <- within_existing %>% rename(test_layer  = layer, in_test  = present)
present_train <- within_existing %>% rename(train_layer = layer, in_train = present)

cross_annot <- df %>%
  filter(train_layer != test_layer) %>%
  left_join(present_test,  by = c("test_layer","node_from","node_to")) %>%
  left_join(present_train, by = c("train_layer","node_from","node_to")) %>%
  mutate(
    # category only defined for test-layer EXISTING links
    category = dplyr::case_when(
      isTRUE(in_test) & isTRUE(in_train) ~ "shared",
      isTRUE(in_test) & !isTRUE(in_train) ~ "unique",
      TRUE ~ NA_character_
    )
  )

# =========================================================
# 3) EVALUATION on removed==1, per (train,test,itr), for unique/shared
#    - binarize predictions with sigmoid + best_optimal_threshold
#    - compute TP/FP/FN and F1
#    We compute F1 for each category by using:
#      positives = removed == 1 & category == {that category}
#      negatives = removed == 1 & original_links == 0
#    (positives from the *other* category are excluded from that category's pool)
# =========================================================
cross_eval <- cross_annot %>%
  filter(removed == 1) %>%
  mutate(
    y_true = as.integer(original_links != 0),
    y_prob = sigmoid(predicted_values),
    y_pred = as.integer(y_prob >= best_optimal_threshold)
  )

metric_fun <- function(d) {
  tp <- sum(d$y_pred == 1 & d$y_true == 1, na.rm = TRUE)
  fp <- sum(d$y_pred == 1 & d$y_true == 0, na.rm = TRUE)
  fn <- sum(d$y_pred == 0 & d$y_true == 1, na.rm = TRUE)
  f1 <- if ((2*tp + fp + fn) > 0) 2*tp/(2*tp + fp + fn) else NA_real_
  tibble(tp = tp, fp = fp, fn = fn, f1 = f1)
}

# F1 for SHARED: use all removed negatives + removed positives tagged "shared"
res_shared <- cross_eval %>%
  filter(y_true == 0 | category == "shared") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(category = "shared")

# F1 for UNIQUE: use all removed negatives + removed positives tagged "unique"
res_unique <- cross_eval %>%
  filter(y_true == 0 | category == "unique") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(category = "unique")

f1_by_itr <- bind_rows(res_shared, res_unique)

# (optional) sanity: how many positives per group?
pos_counts <- cross_eval %>%
  filter(y_true == 1, !is.na(category)) %>%
  count(train_layer, test_layer, itr, category, name = "n_pos_in_cat")

# =========================================================
# 4) Compare F1 across train/test combos for shared vs unique
# =========================================================
# Per (train,test,itr,category)
f1_by_itr  # <- main per-iteration results

# Summary per (train,test,category)
f1_combo_summary <- f1_by_itr %>%
  group_by(train_layer, test_layer, category) %>%
  summarise(
    n_itrs   = n(),
    mean_f1  = mean(f1, na.rm = TRUE),
    median_f1= median(f1, na.rm = TRUE),
    .groups = "drop"
  )

# Difference (shared - unique) per (train,test,itr)
f1_diff_per_itr <- f1_by_itr %>%
  select(train_layer, test_layer, itr, category, f1) %>%
  tidyr::pivot_wider(names_from = category, values_from = f1) %>%
  mutate(diff_shared_minus_unique = shared - unique)

# Difference (shared - unique) averaged per (train,test)
f1_diff_summary <- f1_diff_per_itr %>%
  group_by(train_layer, test_layer) %>%
  summarise(
    mean_diff = mean(diff_shared_minus_unique, na.rm = TRUE),
    median_diff = median(diff_shared_minus_unique, na.rm = TRUE),
    .groups = "drop"
  )

# Useful objects to inspect:
# within_existing, cross_annot (with categories), f1_by_itr, f1_combo_summary, f1_diff_per_itr, f1_diff_summary, pos_counts
