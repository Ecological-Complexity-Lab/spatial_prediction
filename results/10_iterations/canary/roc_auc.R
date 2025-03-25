df <- read_csv('working_df_all_itr_25_binary.csv')


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

# Plot the PR curve with the AUC (integral) in the title
plot(pr_obj, main = paste("PR Curve (AUC =", round(pr_obj$auc.integral, 2), ")"))

library(dplyr)

# Compute ROC AUC for each group (train_layer, test_layer, itr)
roc_results <- df %>%
  group_by(train_layer, test_layer, itr) %>%
  summarize(
    auc_roc = auc(roc(response = original_links, predictor = predicted_prob_sigm))
  )

print(roc_results)
