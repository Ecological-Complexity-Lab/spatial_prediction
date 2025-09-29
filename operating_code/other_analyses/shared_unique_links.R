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
      in_test & in_train ~ "shared",
      in_test & is.na(in_train) ~ "unique_in_P",
      is.na(in_test) & in_train ~ "unique_in_A",
      is.na(in_test) & is.na(in_train) ~ "non_link",
      TRUE ~ NA_character_
    )
  )

table(cross_annot$category)
sum(is.na(cross_annot$category))

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

# Shows counts of 0/1 by category; negatives will appear under NA (expected)
with(cross_eval, table(category, y_true, useNA = "ifany"))

eff_counts_by_group <- cross_eval %>%
  group_by(train_layer, test_layer, itr) %>%
  summarise(
    n_neg        = sum(y_true == 0, na.rm = TRUE),
    n_pos_shared = sum(y_true == 1 & category == "shared", na.rm = TRUE),
    n_pos_unique = sum(y_true == 1 & category == "unique_in_P", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(n_pos_shared, n_pos_unique),
               names_to = "category", values_to = "n_pos") %>%
  mutate(
    category = ifelse(category == "n_pos_shared", "shared", "unique_in_P"),
    n_total  = n_neg + n_pos,
    pos_rate = n_pos / pmax(n_total, 1)
  )

eff_counts_by_group
with(eff_counts_by_group, table(category, n_pos))

metric_fun <- function(d) {
  tp <- sum(d$y_pred == 1 & d$y_true == 1, na.rm = TRUE)
  fp <- sum(d$y_pred == 1 & d$y_true == 0, na.rm = TRUE)
  fn <- sum(d$y_pred == 0 & d$y_true == 1, na.rm = TRUE)
  recall <- (tp / (tp + fn))
  precision <- tp / (tp + fp)
  f1 <- 2 * (precision * recall) / (precision + recall)
  tibble(tp = tp, fp = fp, fn = fn, f1 = f1)
}

# add a step to balance f1 somehow

# F1 for SHARED: use all removed negatives + removed positives tagged "shared"
res_shared <- cross_eval %>%
  filter(y_true == 0 | category == "shared") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(category = "shared")

# F1 for UNIQUE: use all removed negatives + removed positives tagged "unique"
res_unique <- cross_eval %>%
  filter(y_true == 0 | category == "unique_in_P") %>%
  group_by(train_layer, test_layer, itr) %>%
  do(metric_fun(.)) %>%
  ungroup() %>%
  mutate(category = "unique_in_P")

f1_by_itr <- bind_rows(res_shared, res_unique)



f1_by_itr %>%
  mutate(l=tp+fp) %>% 
  group_by(category) %>% 
  summarise(m=mean(tp),
            n=n(),
            s=sum(l)) # we have much mire links considered for unique interactions than shared. so the higher f1 for unique interactions might be size effect.

# ---- compare ----
df_plot_f1 <- f1_by_itr %>%
  filter(category %in% c("unique_in_P","shared"), !is.na(f1)) %>%
  mutate(category = factor(category, levels = c("unique_in_P","shared")))

# --- Basic boxplot across all iters/combos ---
ggplot(df_plot_f1, aes(x = category, y = f1, fill = category)) +
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
    values = c("unique_in_P" = "lightsteelblue2",
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
