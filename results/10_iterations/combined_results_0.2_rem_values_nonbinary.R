# ---- exploring the results df ----
# no aggregation yet

## ---- results so far ----
# binarization based on thresholding (0.5 == 0 otherwise 1) is better than sigmoid. 7.42 trues:falses for thresholding, 0.7992887 trues:falses for sigmoid (we get more falses than trues). probably did it wrong.
# binarization of the matrices *before* imputation reduces the range of the predicted values by approximately 10 times
## ---- to do ----
# find k and lambda for which the results are most accurate
# find the best way to binarize the predicted values
# average the predicted values for the same interactions 
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load results ----
d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv')

## ---- working on one iteration to test ----

# d %>%
#   filter(itr==1) %>%
#   filter(k==2) %>%
#   filter(lambda==0.01) %>%
#   write_csv('test_df_itr_1_25_binary.csv') # create a small test dataframe with 1 iteration and the smallest k and lambda.

# analysing binarization in the test dataframe
d <- read_csv('test_df_itr_1_25_binary.csv')
d %>% 
  ggplot(aes(original_links,predicted_values))+geom_point()

d <- d %>% 
  filter(removed == 1) %>%
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) 

d %>% 
  ggplot(aes(original_links,predicted_prob_sigm))+geom_point()

## ---- some evaluators ----


# Calculate metrics for each unique combination of train_layer and test_layer

specificity = TN / (TN + FP)
precision = TP / (TP + FP)
recall = TP / (TP + FN) # split to recall of removed and non-removed links
f1_score = 2 * (precision * recall) / (precision + recall)
balanced_accuracy = (recall + specificity) / 2
mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

result <- d %>%
  group_by(emln_id, train_layer, test_layer, lambda, k) %>%
  summarise(
    TP = sum(original_links == 1 & predicted_bin_sigm == 1),
    FN = sum(original_links == 1 & predicted_bin_sigm == 0),
    TN = sum(original_links == 0 & predicted_bin_sigm == 0),
    FP = sum(original_links == 0 & predicted_bin_sigm == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  ungroup()

# View the result
print(result)

## ---- visualization ----
heatmap_plot_rem <- ggplot(result, aes(x = train_layer, y = test_layer, fill = recall)) +
  geom_tile(color = "white", linewidth = 0.1) +  # Add white borders for better tile distinction
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  # Set NA values to gray
  labs(title = "recall: with 20% link removal", x = "Training layer", y = "Predicted layer", fill = "Recall") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional)
  ) +
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure each x tick has a label
  scale_y_continuous(breaks = seq(1, 7, by = 1))

print(heatmap_plot_rem)

heatmap_plot_rem <- ggplot(result, aes(x = train_layer, y = test_layer, fill = f1_score)) +
  geom_tile(color = "white", linewidth = 0.1) +  # Add white borders for better tile distinction
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  # Set NA values to gray
  labs(title = "f1_score: with 20% link removal", x = "Training layer", y = "Predicted layer", fill = "f1_score") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank()   # Remove panel background (optional)
  ) +
  scale_x_continuous(breaks = seq(1, 7, by = 1)) +  # Ensure each x tick has a label
  scale_y_continuous(breaks = seq(1, 7, by = 1))

print(heatmap_plot_rem)
