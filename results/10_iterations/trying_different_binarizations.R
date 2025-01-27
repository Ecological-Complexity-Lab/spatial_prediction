# ---- exploring the results df ----
# no aggregation yet

## ---- results so far ----
# binarization based on thresholding (0.5 == 0 otherwise 1) is better than sigmoid. 7.42 trues:falses for thresholding, 0.7992887 trues:falses for sigmoid (we get more falses than trues). probably did it wrong.
# binarization of the matrices *before* imputation reduces the range of the predicted values by approximately 10 times
# the variation in predicting zeros is consistantly smaller than predicting ones
## ---- to do ----
# find k and lambda for which the results are most accurate
# find the best way to binarize the predicted values
# average the predicted values for the same interactions 
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)
library(emln)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Min-max normalization function
normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply normalization to predicted_value
d <- d %>%
  mutate(probability = normalize_min_max(predicted_value))


## ---- load results ----
d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv')

## ---- working on one iteration to test ----

d %>%
  filter(itr==1) %>%
  filter(k==2) %>%
  filter(lambda==0.01) %>%
  write_csv('test_df_itr_1_60_binary.csv') # create a small test dataframe with 1 iteration and the smallest k and lambda.

# analysing binarization in the test dataframe
d <- read_csv('test_df_itr_1_25_binary.csv')
d %>% 
  ggplot(aes(original_links,predicted_values))+geom_point()

# alternatively:
d_bin_1 <- d %>% 
  filter(removed == 1) %>%
  #mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_1 = if_else(predicted_values > 0, 1, 0)) 

# or use min-max normalization:
d_minmax <- d %>% 
  filter(removed == 1) %>%
  mutate(predicted_prob_min_max = normalize_min_max(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_minmax = if_else(predicted_values > 0.5, 1, 0))

## ---- some evaluators ----

# Calculate metrics for each unique combination of train_layer and test_layer

result_bin_1 <- d_bin_1 %>%
  group_by(emln_id, train_layer, test_layer, lambda, k) %>%
  summarise(
    TP = sum(original_links == 1 & predicted_bin_1 == 1),
    FN = sum(original_links == 1 & predicted_bin_1 == 0),
    TN = sum(original_links == 0 & predicted_bin_1 == 0),
    FP = sum(original_links == 0 & predicted_bin_1 == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  ungroup()

result_minamax <- d_minmax %>%
  group_by(emln_id, train_layer, test_layer, lambda, k) %>%
  summarise(
    TP = sum(original_links == 1 & predicted_bin_1 == 1),
    FN = sum(original_links == 1 & predicted_bin_1 == 0),
    TN = sum(original_links == 0 & predicted_bin_1 == 0),
    FP = sum(original_links == 0 & predicted_bin_1 == 1),
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

# Get layer names
net <- emln::load_emln(25)
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Display the updated tibble
print(net_name)

# Perform left joins to replace IDs with names
result <- result %>%
  mutate(train_layer = as.integer(as.character(train_layer)),  # Ensure train_layer is an integer
         test_layer = as.integer(as.character(test_layer))) %>% # Ensure test_layer is an integer
  left_join(net_name, by = c("train_layer" = "layer_id")) %>% 
  rename(train_layer_name = name) %>%                          # Use a different name for clarity
  left_join(net_name, by = c("test_layer" = "layer_id")) %>%  
  rename(test_layer_name = name)                               # Use a different name for clarity
# Rename joined column


# Set factor levels for training and predicting layers
layer_levels <- as.character(1:7)
result$train_layer <- factor(result$train_layer, levels = layer_levels)
result$test_layer <- factor(result$test_layer, levels = unique(result$test_layer))


result$diagonal <- result$train_layer == result$test_layer
layer_to_layer_plot <- 
  ggplot(result, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result[result$train_layer == result$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "f1 score") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

print(layer_to_layer_plot)

pdf(paste0(output_folder,'pr_layer_to_layer.pdf'), 7, 7)
print(layer_to_layer_plot)
dev.off()

layer_to_layer_plot_canary <- 
  ggplot(result, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result[result$train_layer == result$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "f1 score") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

print(layer_to_layer_plot_canary)

# check the effect of k and lambda on a test tibble with 1 iteration

d <- read_csv('nonbinary_equal_0_1_removal_60_1.csv')

## ---- moving to all iterations ----
