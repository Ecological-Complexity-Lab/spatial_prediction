### ---- selecting optimal threshold ----
library(dplyr)
library(tidyr)
library(ggplot2)

# 0) your thresholds
thresholds <- seq(0, 1, by = 0.1)

# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

# 2) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(emln_id, train_layer, test_layer, itr, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy= (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) /
      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)#,
    #.groups = "drop"
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup() 

# 3) average across your emln_id/layer combos and pivot long
df_avg <- df_thresh %>%
  group_by(threshold) %>%
  summarise(across(
    c(specificity, precision, recall,
      f1_score, balanced_accuracy, mcc),
    mean, na.rm = TRUE
  )) %>%
  pivot_longer(-threshold,
               names_to  = "metric",
               values_to = "value")

# 4) plot
ggplot(df_avg, aes(threshold, value, color = metric)) +
  geom_line(size = 1) +
  labs(
    x     = "Probability threshold",
    y     = "Average metric",
    color = "Metric"
  ) +
  tme

library(dplyr)
library(tidyr)

# 1) pivot to wide so F1 and balanced_accuracy are columns
df_wide <- df_avg %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  arrange(threshold)

# 2) discrete approx: minimize abs difference
best_discrete <- df_wide %>%
  mutate(absdiff = abs(f1_score - balanced_accuracy)) %>%
  slice_min(absdiff, n = 1)

# results:
best_discrete_threshold <- best_discrete$threshold
best_discrete_f1        <- best_discrete$f1_score
best_discrete_bal_acc   <- best_discrete$balanced_accuracy

# 3) linear interpolation between the two points that straddle the crossing
# compute difference vector
diff_vec <- df_wide$f1_score - df_wide$balanced_accuracy

# find the first place where it changes sign
ix <- which(diff(sign(diff_vec)) != 0)[1]

if (!is.na(ix)) {
  # the two bounding thresholds
  t1 <- df_wide$threshold[ix]
  t2 <- df_wide$threshold[ix + 1]
  # their metric values
  f1_1 <- df_wide$f1_score[ix];  f1_2 <- df_wide$f1_score[ix + 1]
  ba_1 <- df_wide$balanced_accuracy[ix];  ba_2 <- df_wide$balanced_accuracy[ix + 1]
  # slopes
  s1 <- (f1_2 - f1_1) / (t2 - t1)
  s2 <- (ba_2 - ba_1) / (t2 - t1)
  # solve for t where f1_1 + s1*(t - t1) = ba_1 + s2*(t - t1)
  t_intersect <- t1 + (ba_1 - f1_1) / (s1 - s2)
} else {
  t_intersect <- NA_real_
}

# now you have:
best_discrete_threshold  #: the threshold with minimal |F1 - bal_acc|
t_intersect              #: the linearly‐interpolated exact crossing
