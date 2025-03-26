df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

df <- df %>%
  filter(removed == 1) %>% mutate(predicted_prob_sigm = sigmoid(predicted_values))

robust_sigmoid <- function(x, center = median(x), scale = mad(x)) {
  1 / (1 + exp(-(x - center) / scale))
}

df <- df %>%
 filter(removed == 1) %>% mutate(predicted_prob_sigm = robust_sigmoid(predicted_values))

library(pROC)
library(PRROC)

# ROC curve
roc_obj <- roc(response = df$original_links, predictor = df$predicted_prob_sigm)
auc_roc <- auc(roc_obj)

# Plot the ROC curve with the AUC in the title
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_roc, 2), ")"))

# pr curve
# For the PR curve, we separate the scores for the positive class (original_links==1) and the negative class (original_links==0)
pr_obj <- pr.curve(
  scores.class0 = df$predicted_prob_sigm[df$original_links == 1],
  scores.class1 = df$predicted_prob_sigm[df$original_links == 0],
  curve = TRUE
)

# Calculate the positive class ratio
pos_ratio <- sum(df$original_links == 1) / nrow(df)

# Plot the PR curve with the AUC (integral) in the title
plot(pr_obj, main = paste("PR Curve (AUC =", round(pr_obj$auc.integral, 2), ")"))

# Add the random guess line (horizontal line at the positive ratio)
abline(h = pos_ratio, col = "red", lty = 2)


df %>% group_by(original_links) %>% summarise(n=n()) %>% arrange(desc(n))


# Additional analysis to better understand model performance
library(ggplot2)

# Examine distribution of predictions by true class
ggplot(df, aes(x = predicted_prob_sigm, fill = factor(original_links))) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Predicted Probabilities by True Class",
       x = "Predicted Probability", y = "Density",
       fill = "Original Link") +
  theme_minimal()












# 1 itr

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


