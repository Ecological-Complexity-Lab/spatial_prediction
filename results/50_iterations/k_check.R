# ---- check performance according to k ----
## ---- functions ----
# transformations

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- data ----
df_50_site <- read_csv('canary_weighted_scaled_site_net_60_50itr.csv')

## ---- evaluation ----
df_removed <- df_50_site %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))


# result_summary <- df_removed %>%
#   group_by(emln_id, train_layer, test_layer, k, itr) %>%
#   summarise(
#     TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
#     FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
#     TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
#     FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
#     specificity = TN / (TN + FP),
#     precision = TP / (TP + FP),
#     recall = TP / (TP + FN),
#     f1_score = 2 * (precision * recall) / (precision + recall),
#     balanced_accuracy = (recall + specificity) / 2,
#     mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
#     mse = mean((predicted_values - original_links)^2, na.rm = TRUE),
#     rmse = sqrt(mse)
#   ) %>%
#   ungroup() %>%
#   group_by(emln_id, train_layer, test_layer, k) %>%
#   summarise(
#     TP = mean(TP, na.rm = TRUE),
#     FN = mean(FN, na.rm = TRUE),
#     TN = mean(TN, na.rm = TRUE),
#     FP = mean(FP, na.rm = TRUE),
#     specificity = mean(specificity, na.rm = TRUE),
#     precision = mean(precision, na.rm = TRUE),
#     recall = mean(recall, na.rm = TRUE),
#     f1_score = mean(f1_score, na.rm = TRUE),
#     balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
#     mcc = mean(mcc, na.rm = TRUE),
#     mse = mean(mse, na.rm = TRUE),
#     rmse = mean(rmse, na.rm = TRUE)
#   ) %>%
#   ungroup()

print(head(df_50_site, 55), width = Inf)

result_summary <- df_removed %>%
  group_by(emln_id, train_layer, test_layer, k, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin_sigm == 1),
    FN = sum(original_binary == 1 & predicted_bin_sigm == 0),
    TN = sum(original_binary == 0 & predicted_bin_sigm == 0),
    FP = sum(original_binary == 0 & predicted_bin_sigm == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
    mse = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)
  ) %>%
  ungroup()

# the results for all ks are the same...
