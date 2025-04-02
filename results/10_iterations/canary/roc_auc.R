# ---- libraries ----
library(pROC)
library(PRROC)
library(tidyverse)
library(ggplot2)

# ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))

# ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

robust_sigmoid <- function(x, center = median(x), scale = mad(x)) {
  1 / (1 + exp(-(x - center) / scale))
}

normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

tanh_transform <- function(x){
  (tanh(x)+1)/2
}

clip_transform <- function(x){
  case_when(x<0~0,
            x>1~1,
            TRUE~x)
}
# ---- data ----
df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv') # binary version, filtered k and lambda
df <- read.csv('aggregated_equal_0_1_removal_60_1.csv') # binary version, no filtering
df <- read.csv('nonbinary_equal_0_1_removal_60_0.csv') # weighted version
df
df <- read.csv('binary_equal_0_1_removal_scaling_60_1.csv') # scaled
# df <- read.csv('binary_equal_0_1_removal_scaling_lambda_5_60_1.csv') # test
# ---- analysis ----

# lambdas=unique(df$lambda)
#0.01
# df <- df %>% filter(k == 2 & lambda == lambdas[1]) # if we want specific k and lambda
# df <- df %>% filter(k == 2)
# df <- df %>% filter(itr == 1) # if we want to look at 1 zero removal iteration

df <- df %>%
  filter(removed == 1) %>%
  mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>% 
  mutate(predicted_prob_robust_sigm = robust_sigmoid(predicted_values)) %>% 
  mutate(predicted_minmax = normalize_min_max(predicted_values)) %>% 
  mutate(predicted_tanh = tanh_transform(predicted_values)) %>% 
  mutate(predicted_clip = clip_transform(predicted_values))

# # ---- roc pr ----
# # ROC curve
# roc_obj <- roc(response = df$original_links, predictor = df$predicted_prob_sigm)
# auc_roc <- auc(roc_obj)
# 
# # Plot the ROC curve with the AUC in the title
# plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_roc, 2), ")"))
# 
# # pr curve
# # For the PR curve, we separate the scores for the positive class (original_links==1) and the negative class (original_links==0)
# pr_obj <- pr.curve(
#   scores.class0 = df$predicted_prob_sigm[df$original_links == 1],
#   scores.class1 = df$predicted_prob_sigm[df$original_links == 0],
#   curve = TRUE
# )
# 
# # Calculate the positive class ratio
# pos_ratio <- sum(df$original_links == 1) / nrow(df)
# 
# # Plot the PR curve with the AUC (integral) in the title
# plot(pr_obj, main = paste("PR Curve (AUC =", round(pr_obj$auc.integral, 2), ")"))
# 
# # Add the random guess line (horizontal line at the positive ratio)
# abline(h = pos_ratio, col = "red", lty = 2)
# 
# 
# df %>% group_by(original_links) %>% summarise(n=n()) %>% arrange(desc(n))
# 
# # ---- test different transformations ----

# Function to plot ROC curve with ggplot2
plot_roc_curve <- function(true_labels, predicted_scores) {
  # Create the ROC object and compute AUC
  roc_obj <- roc(response = true_labels, predictor = predicted_scores)
  auc_val <- auc(roc_obj)
  
  # Build a data frame from the ROC object for ggplot2
  df_roc <- data.frame(
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
  
  # Generate the ROC plot
  p <- ggplot(df_roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # ROC curve line
    geom_abline(intercept = 0, slope = 1,                       # Diagonal line (random classifier)
                linetype = "dashed", color = "salmon") +
    labs(title = paste("ROC curve (AUC =", round(auc_val, 2), ")"),
         x = "False positive rate", y = "True positive rate") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

# Function to plot PR curve with ggplot2
plot_pr_curve <- function(true_labels, predicted_scores) {
  # Separate scores by class: positive (true label==1) and negative (true label==0)
  scores_pos <- predicted_scores[true_labels == 1]
  scores_neg <- predicted_scores[true_labels == 0]
  
  # Create the PR curve object; curve=TRUE returns the full curve data
  pr_obj <- pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = TRUE)
  
  # Calculate the positive class ratio for the random baseline line
  pos_ratio <- sum(true_labels == 1) / length(true_labels)
  
  # Convert the curve matrix to a data frame and set column names
  df_pr <- as.data.frame(pr_obj$curve)
  colnames(df_pr) <- c("recall", "precision", "threshold")
  
  # Generate the PR plot
  p <- ggplot(df_pr, aes(x = recall, y = precision)) +
    geom_line(color = "lightsteelblue", size = 1) +                      # PR curve line
    geom_hline(yintercept = pos_ratio,                          # Baseline: random classifier performance
               linetype = "dashed", color = "salmon") +
    labs(title = paste("PR curve (AUC =", round(pr_obj$auc.integral, 2), ")"),
         x = "Recall", y = "Precision") +
    theme_minimal() + tme                                           # Clean theme
  print(p)
}

# mutate(original_links_binary = ifelse(original_links == 0, 0, 1)) %>% 
#   mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>% 
#   mutate(predicted_prob_robust_sigm = robust_sigmoid(predicted_values)) %>% 
#   mutate(predicted_minmax = normalize_min_max(predicted_values)) %>% 
#   mutate(predicted_tanh = tanh_transform(predicted_values)) %>% 
#   mutate(predicted_clip = clip_transform(predicted_values))

# For the ROC curve:
plot_roc_curve(df$original_links, df$predicted_values)

# For the PR curve:
plot_pr_curve(df$original_links, df$predicted_values)


# ---- distribution of predictions by true class ---- 
ggplot(df, aes(x = predicted_prob_sigm, fill = factor(original_links))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of predicted probabilities \nby true class",
       x = "Predicted probability", y = "Density",
       fill = "Original Link") +
  theme_minimal() + tme

# for weights
ggplot(df, aes(x = original_links, y = predicted_prob_sigm)) +
  geom_point(alpha = 0.6, size = 2, color = "lightsteelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Observed weight",
    y = "Predicted probability"
  ) +
  theme_minimal() + tme

# ---- 1 itr ----

df_itr <- df %>%
  filter(itr == 1)

# ROC curve
roc_obj_itr <- roc(response = df_itr$original_links, predictor = df_itr$predicted_prob_sigm)
auc_roc_itr <- auc(roc_obj_itr)

# Plot the ROC curve with the AUC in the title
plot(roc_obj_itr, main = paste("ROC Curve (AUC =", round(auc_roc, 2), ")"))

# pr curve
# For the PR curve, we separate the scores for the positive class (original_links==1) and the negative class (original_links==0)
pr_obj_itr <- pr.curve(
  scores.class0 = df_itr$predicted_prob_sigm[df_itr$original_links == 1],
  scores.class1 = df_itr$predicted_prob_sigm[df_itr$original_links == 0],
  curve = TRUE
)

# Plot the PR curve with the AUC (integral) in the title
plot(pr_obj_itr, main = paste("PR Curve (AUC =", round(pr_obj$auc.integral, 2), ")"))


best_threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
print(best_threshold)

library(dplyr)

# Compute ROC AUC for each group (train_layer, test_layer, itr)
roc_results <- df %>%
  group_by(train_layer, test_layer, itr) %>%
  summarize(
    auc_roc = auc(roc(response = original_links, predictor = predicted_prob_sigm))
  )

print(roc_results)

thresholds <- seq(0, 1, by = 0.01)
f1_scores <- sapply(thresholds, function(th) {
  predicted <- ifelse(df$predicted_prob_sigm >= th, 1, 0)
  precision <- sum(predicted == 1 & df$original_links == 1) / sum(predicted == 1)
  recall <- sum(predicted == 1 & df$original_links == 1) / sum(df$original_links == 1)
  if (is.na(precision) || is.na(recall) || (precision + recall) == 0) return(0)
  2 * precision * recall / (precision + recall)
})
best_threshold_pr <- thresholds[which.max(f1_scores)]
print(best_threshold_pr)


