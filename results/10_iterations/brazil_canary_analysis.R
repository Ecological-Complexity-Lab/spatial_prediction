# ---- exploring the results df ----
# no aggregation yet

## ---- results so far ----
# binarization based on thresholding (0.5 == 0 otherwise 1) is better than sigmoid. 7.42 trues:falses for thresholding, 0.7992887 trues:falses for sigmoid (we get more falses than trues). probably did it wrong.
# binarization of the matrices *before* imputation reduces the range of the predicted values by approximately 10 times
# the variation in predicting zeros is consistantly smaller than predicting ones
# k=2 is much better for network 25. lambda values result in similar f1 values but lambda 0.1 is a bit better. also for the canary islands.
# results for averaging evaluators for all iterations are suspiciously similar to the evaluators for one iteration...
# but the results for different iterations are not the same
## ---- to do ----
# find the best way to binarize the predicted values
# average the predicted values for the same interactions 
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)
library(emln)
library(reshape2)
library(ggpubr)

tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load results ----
d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv')

## ---- working on one iteration to test ----

d %>%
  filter(itr==1) %>%
  filter(k==2) %>%
  filter(lambda==0.01) %>%
  write_csv('test_df_itr_1_60_binary.csv') # create a small test dataframe with 1 iteration and the smallest k and lambda.

# analysing binarization in the test dataframe
d <- read_csv('test_df_itr_1_60_binary.csv')
d %>% 
  ggplot(aes(original_links,predicted_values))+geom_point()

d <- d %>% 
  filter(removed == 1) %>%
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) %>% 
  mutate(predicted_bin_over0 = ifelse(predicted_values > 0, 1, 0))

d %>% 
  ggplot(aes(original_links,predicted_prob_sigm))+geom_point()

## ---- some evaluators for 1 iteration ----


# Calculate metrics for each unique combination of train_layer and test_layer

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

## ---- comparing binarization using sigmoid and above 0 ----

result <- d %>%
  # Only keep rows that have removed == 1 if desired
  # filter(removed == 1) %>%
  group_by(emln_id, train_layer, test_layer, lambda, k) %>%
  summarise(
    # --- For predicted_bin_sigm ---
    TP_sigm = sum(original_links == 1 & predicted_bin_sigm == 1),
    TN_sigm = sum(original_links == 0 & predicted_bin_sigm == 0),
    FP_sigm = sum(original_links == 0 & predicted_bin_sigm == 1),
    FN_sigm = sum(original_links == 1 & predicted_bin_sigm == 0),
    
    precision_sigm = TP_sigm / (TP_sigm + FP_sigm),
    recall_sigm    = TP_sigm / (TP_sigm + FN_sigm),
    f1_sigm = 2 * (precision_sigm * recall_sigm) / (precision_sigm + recall_sigm),
    
    # Optional: specificity, balanced_accuracy, etc. as in your example
    specificity_sigm = TN_sigm / (TN_sigm + FP_sigm),
    balanced_acc_sigm = (recall_sigm + specificity_sigm) / 2,
    
    # --- For predicted_bin_over0 ---
    TP_over0 = sum(original_links == 1 & predicted_bin_over0 == 1),
    TN_over0 = sum(original_links == 0 & predicted_bin_over0 == 0),
    FP_over0 = sum(original_links == 0 & predicted_bin_over0 == 1),
    FN_over0 = sum(original_links == 1 & predicted_bin_over0 == 0),
    
    precision_over0 = TP_over0 / (TP_over0 + FP_over0),
    recall_over0    = TP_over0 / (TP_over0 + FN_over0),
    f1_over0 = 2 * (precision_over0 * recall_over0) / (precision_over0 + recall_over0),
    
    specificity_over0 = TN_over0 / (TN_over0 + FP_over0),
    balanced_acc_over0 = (recall_over0 + specificity_over0) / 2,
    
    .groups = "drop"   # Stop grouping after summarise
  )

f1_summary <- result %>%
  summarise(
    max_f1_sigm = max(f1_sigm, na.rm = TRUE),
    avg_f1_sigm = mean(f1_sigm, na.rm = TRUE),
    
    max_f1_over0 = max(f1_over0, na.rm = TRUE),
    avg_f1_over0 = mean(f1_over0, na.rm = TRUE)
  )

f1_summary


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
layer_levels <- as.character(1:7) # for Brazil
layer_levels <- as.character(1:14) # for Canary Islands
result$train_layer <- factor(result$train_layer, levels = layer_levels)
result$test_layer <- factor(result$test_layer, levels = unique(result$test_layer))


result$diagonal <- result$train_layer == result$test_layer
layer_to_layer_plot <- 
  ggplot(result, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result[result$train_layer == result$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "balanced accuracy") +
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
  ggplot(result, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result[result$train_layer == result$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "balanced accuracy") +
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

d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv') # brazil
d <- read_csv('nonbinary_equal_0_1_removal_60_1.csv') # canary islands

d <- d %>%
  filter(removed == 1) %>% 
  filter(k == 2) %>% 
  filter(lambda == 0.1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) #%>%
  #write_csv('working_df_all_itr_25_binary.csv')

#d <- read_csv('working_df_all_itr_25_binary.csv')

result_summary <- d %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
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
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  ungroup()

# Get layer names
net <- emln::load_emln(25) # brazil
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Display the updated tibble
print(net_name)

# Perform left joins to replace IDs with names
result_summary <- result_summary %>%
  mutate(train_layer = as.integer(as.character(train_layer)),  # Ensure train_layer is an integer
         test_layer = as.integer(as.character(test_layer))) %>% # Ensure test_layer is an integer
  left_join(net_name, by = c("train_layer" = "layer_id")) %>% 
  rename(train_layer_name = name) %>%                          # Use a different name for clarity
  left_join(net_name, by = c("test_layer" = "layer_id")) %>%  
  rename(test_layer_name = name)                               # Use a different name for clarity

#result_summary <- result_summary %>% write_csv('summary_df_all_itr_1_25_binary.csv')

# Set factor levels for training and predicting layers
#layer_levels <- as.character(1:7) # for Brazil
# layer_levels <- as.character(1:14) # for Canary Islands
# result_summary$train_layer <- factor(result_summary$train_layer, levels = layer_levels)
# result_summary$test_layer <- factor(result_summary$test_layer, levels = unique(result_summary$test_layer))


result_summary$diagonal <- result_summary$train_layer == result_summary$test_layer

#write_csv(result_summary, 'working_df_all_itr_60_binary_names.csv')

layer_to_layer_plot_all_itr_brzail <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "f1") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

print(layer_to_layer_plot_all_itr_brzail)

# making a half heatmap

result_summary_filtered <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

layer_to_layer_plot_all_itr_brzail_half <-  
  ggplot(result_summary_filtered, aes(x = train_layer, y = test_layer, fill = f1_score)) +
  geom_tile(color = "white", linewidth = 0.1) +  
  geom_tile(data = result_summary_filtered[result_summary_filtered$train_layer == result_summary_filtered$test_layer, ],
            color = "black", linewidth = 1.2) +  
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +  
  scale_x_discrete(limits = rev(unique(result_summary$train_layer))) +  # Reverse x-axis order
  scale_y_discrete(limits = rev(unique(result_summary$test_layer))) +
  labs(x = "Training layer", y = "Predicted layer", fill = "f1") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  
  ) +
  coord_fixed()
layer_to_layer_plot_all_itr_brzail_half


result_summary_filtered <- read.csv('working_df_all_itr_25_binary_names_half_matrix.csv')

# Now you should have two factor columns that won't trigger duplicate-level errors,
# assuming each numeric layer truly lines up with exactly one layer name.
result_summary_filtered <- result_summary_filtered %>%
  mutate(
    # Make them factors, reversing the sorted numeric vector
    train_layer_factor = factor(
      train_layer,
      levels = rev(sort(unique(train_layer)))  # largest -> smallest
    ),
    test_layer_factor = factor(
      test_layer,
      levels = rev(sort(unique(test_layer)))
    )
  )

result_summary_filtered <- result_summary_filtered %>%
  mutate(
    train_layer = factor(train_layer, levels = rev(sort(unique(train_layer)))),
    test_layer  = factor(test_layer,  levels = rev(sort(unique(test_layer))))
  )

layer_to_layer_plot_all_itr_brazil_half <- ggplot(
  result_summary_filtered,
  aes(x = train_layer, y = test_layer
      , fill = f1_score)
) +
  geom_tile(color = "white", linewidth = 0.1) +
  
  # If you want black borders on the diagonal:
  geom_tile(
    data = subset(result_summary_filtered, train_layer == test_layer),
    color = "black", linewidth = 1.2
  ) +
  #scale_x_discrete(limits = levels(rev(result_summary_filtered$train_layer_factor)), drop = FALSE) +
  #scale_y_discrete(limits = levels(result_summary_filtered$test_layer_factor), drop = FALSE) +
  #scale_y_continuous(trans = "reverse", breaks = seq(1, by = 1)) +

  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +
  labs(x = "Training layer", y = "Predicted layer", fill = "f1") +
  theme_minimal() +
  theme(
    plot.margin      = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed()

layer_to_layer_plot_all_itr_brazil_half


layer_to_layer_plot_all_itr_brzail_half <- ggplot(
  result_summary_filtered, 
  aes(x = train_layer_name, y = test_layer_name, fill = f1_score)
) +
  geom_tile(color = "white", linewidth = 0.1) +
  # Black borders only for diagonal tiles
  geom_tile(
    data = result_summary_filtered[result_summary_filtered$train_layer == result_summary_filtered$test_layer, ],
    color = "black", linewidth = 1.2
  ) +
  scale_fill_gradient(low = "skyblue", high = "orchid4", na.value = "gray") +
  labs(x = "Training layer", y = "Predicted layer", fill = "f1") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed()

layer_to_layer_plot_all_itr_brzail_half

# canaries 

layer_to_layer_plot_all_itr_canary <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "balanced accuracy") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

print(layer_to_layer_plot_all_itr_canary)

# now make half a heatmap

result_summary_filtered_canary <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

layer_to_layer_plot_all_itr_canary_half <-  
  ggplot(result_summary_filtered_canary, 
         aes(x = factor(train_layer, levels = unique(train_layer)),  
             y = factor(test_layer, levels = rev(unique(test_layer))),  # ✅ FIXED: Closed parentheses here
             fill = f1_score)) +  
  geom_tile(color = "white", linewidth = 0.1) +  
  geom_tile(data = result_summary_filtered_canary[result_summary_filtered_canary$train_layer == result_summary_filtered_canary$test_layer, ], 
            color = "black", linewidth = 1.2) +  
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  
  labs(x = "Training layer", y = "Predicted layer", fill = "f1") +  
  theme_minimal() +  
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  
    panel.background = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  
  ) +  
  coord_fixed()  # ✅ FIXED: Removed extra closing parenthesis


print(layer_to_layer_plot_all_itr_canary_half)

# check if the differences are significant

# checking which k and lambda are the best 
d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv')

d <- d %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) 

result_summary <- d %>%
  group_by(emln_id, train_layer, test_layer, lambda, k, itr) %>%
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
  ungroup() %>%
  group_by(lambda, k) %>%  # Group by lambda and k to average the metrics
  summarise(
    avg_specificity = mean(specificity, na.rm = TRUE),
    avg_precision = mean(precision, na.rm = TRUE),
    avg_recall = mean(recall, na.rm = TRUE),
    avg_f1_score = mean(f1_score, na.rm = TRUE),
    avg_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    avg_mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  ungroup()

# Create a boxplot for average F1 scores across different k values
ggplot(result_summary, aes(x = factor(k), y = avg_f1_score)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(
    title = "Comparison of Average F1 Scores for Different k Values",
    x = "k value",
    y = "Average f1 score"
  ) +
  theme_minimal()

# Create a boxplot for average F1 scores across different lambda values
ggplot(result_summary, aes(x = factor(lambda), y = avg_f1_score)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(
    title = "Comparison of Average F1 Scores for Different lambda Values",
    x = "Lambda value",
    y = "Average f1 score"
  ) +
  theme_minimal()

# checking for the canary islands

d <- read_csv('nonbinary_equal_0_1_removal_60_1.csv')

d <- d %>%
  filter(removed == 1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) 

result_summary <- d %>%
  group_by(emln_id, train_layer, test_layer, lambda, k, itr) %>%
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
  ungroup() %>%
  group_by(lambda, k) %>%  # Group by lambda and k to average the metrics
  summarise(
    avg_specificity = mean(specificity, na.rm = TRUE),
    avg_precision = mean(precision, na.rm = TRUE),
    avg_recall = mean(recall, na.rm = TRUE),
    avg_f1_score = mean(f1_score, na.rm = TRUE),
    avg_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    avg_mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  ungroup()


# Create a boxplot for average F1 scores across different k values
ggplot(result_summary, aes(x = factor(k), y = avg_f1_score)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(
    title = "Comparison of Average F1 Scores for Different k Values",
    x = "k value",
    y = "Average f1 score"
  ) +
  theme_minimal()

# Create a boxplot for average F1 scores across different lambda values
ggplot(result_summary, aes(x = factor(lambda), y = avg_f1_score)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  labs(
    title = "Comparison of Average F1 Scores for Different lambda Values",
    x = "Lambda value",
    y = "Average f1 score"
  ) +
  theme_minimal()

## ---- correlate evaluators with distance ----------
# Load the pairwise distance matrix
distance_matrix <- read.csv("brazil_geographic_distances.csv", row.names = 1)
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- colnames(distance_matrix) <- 1:7

result_summary$train_layer <- as.factor(result_summary$train_layer)
result_summary$test_layer <- as.factor(result_summary$test_layer)

# Create a column with the geographic distance between the layers for each row
result_summary <- result_summary %>%
  mutate(geographic_distance = mapply(function(train, test) distance_matrix[as.character(train), as.character(test)],
                                      train_layer, test_layer))
#result_summary[1:15, c("train_layer","test_layer","geographic_distance")]

correlation <- cor.test(result_summary$balanced_accuracy, result_summary$geographic_distance, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot <- 
  ggplot(result_summary, aes(x = geographic_distance, y = balanced_accuracy)) +
  geom_point(color = "blue", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "purple") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "Balanced accuracy") +
  tme + 
  annotate("text",
           x = 7, y = 0.6,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")

print(cor_plot)

# try f1...
correlation <- cor.test(result_summary$f1_score, result_summary$geographic_distance, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot <- 
  ggplot(result_summary, aes(x = geographic_distance, y = f1_score)) +
  geom_point(color = "blue", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "purple") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "F1 score") +
  tme + 
  annotate("text",
           x = 7, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")

print(cor_plot)

# trying the canary islands...
# result_summary %>% 
#   write_csv('result_summary_canary_with_distance.csv') # summary of all iterations
# 
# # Load the pairwise distance matrix
# distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL) # we need to match the layer names to the numbering in the results table
# 
# # For each row's train_layer ID, find which row in net_name has the same layer_id,
# # and pull out the corresponding 'name'.
# result_summary$train_name <- net_name$name[ match(result_summary$train_layer, net_name$layer_id) ]
# 
# # Same for test_layer
# result_summary$test_name <- net_name$name[ match(result_summary$test_layer, net_name$layer_id) ]
# 
# distance_table_sym <- distance_table %>%
#   bind_rows(
#     distance_table %>%
#       rename(to = from, from = to)  # swap columns
#   ) %>%
#   distinct(from, to, distance_km)
# 
# # result_summary has train_name, test_name, 
# # plus an old "distance_km" you want to replace
# 
# result_summary <- result_summary %>%
#   select(-distance_km) %>%                # remove old distance_km if it exists
#   left_join(
#     distance_table_sym,
#     by = c("train_name" = "from", "test_name" = "to")
#   ) %>%
#   mutate(distance_km = if_else(train_name == test_name,
#                                0,              # distance = 0 if same site
#                                distance_km))   # otherwise, keep joined distance
# 
# 
# 
# result_summary %>%
#   anti_join(distance_table_sym, by = c("train_name" = "from", "test_name" = "to"))

result_summary <- read_csv('result_summary_canary_with_distance.csv')

correlation <- cor.test(result_summary$f1_score, result_summary$distance_km, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, format = "f", digits = 5)

#p_value <- round(correlation$p.value, 5)

label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot_canary <- 
  ggplot(result_summary, aes(x = distance_km, y = f1_score)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "F1 score") +
  annotate("text",
           x = 370, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 3,
           color = "black") +
  tme

print(cor_plot_canary)

# now for half the matrix

result_summary_filtered_canary <- result_summary %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer | train_layer == test_layer)

correlation <- cor.test(result_summary_filtered_canary$f1_score, result_summary_filtered_canary$distance_km, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 5)  # or round as you prefer
p_value <- formatC(correlation$p.value, format = "f", digits = 5)

p_value <- round(correlation$p.value, 3)

label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot_canary_half <- 
  ggplot(result_summary_filtered_canary, aes(x = distance_km, y = f1_score)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "F1 score") +
  annotate("text",
           x = 385, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 3,
           color = "black") +
  tme

print(cor_plot_canary_half)


# try f1...
correlation <- cor.test(result_summary$f1_score, result_summary$geographic_distance, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot <- 
  ggplot(result_summary, aes(x = geographic_distance, y = f1_score)) +
  geom_point(color = "blue", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "purple") +  # Linear regression line
  labs(x = "Geographic distance (km)",
       y = "F1 score") +
  tme + 
  annotate("text",
           x = 7, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")

print(cor_plot)

## ---- check if the predictions of the off-diagonals are different ----

result_table_site <- read.csv('working_df_all_itr_60_binary_names.csv') # site scale
result_table_isl <- read.csv('working_df_island_distance_fidelity_jaccard_netsize.csv') # island scale

result_table_site <- result_table_site %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Diagonal",
    train_layer > test_layer ~ "Off-diagonal (train > test)",
    train_layer < test_layer ~ "Off-diagonal (train < test)"
  ))
result_table_isl <- result_table_isl %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Diagonal",
    train_layer > test_layer ~ "Off-diagonal (train > test)",
    train_layer < test_layer ~ "Off-diagonal (train < test)"
  ))

# Define custom colors for each group
custom_colors <- c("Diagonal" = "#1b9e77", 
                   "Off-diagonal (train > test)" = "#d95f02", 
                   "Off-diagonal (train < test)" = "#7570b3")

# Function to create notched boxplots with customizations
plot_f1_boxplot <- function(data, metric, y_axis_label = "Balanced accuracy") {
  ggplot(data, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = FALSE, alpha = 0.4, color = "black") +  # Notched, semi-transparent, black outline
    theme_minimal() +
    labs(x = "Layer comparison",
         y = y_axis_label) +  # Custom y-axis title
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black frame around the plot
    ) +
    tme +
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    stat_compare_means(method = "t.test", label = "p.signif", hide.ns = FALSE, comparisons = list(
      c("Diagonal", "Off-diagonal (train > test)"),
      c("Diagonal", "Off-diagonal (train < test)"),
      c("Off-diagonal (train > test)", "Off-diagonal (train < test)")
    ))  # Pairwise t-tests with significance labels
}

# Generate boxplots for each F1 score metric
island_f1 <- plot_f1_boxplot(result_table_isl, metric = "balanced_accuracy")

# combined scales plot
site_f1 <- plot_f1_boxplot(result_table_site, metric = "balanced_accuracy")

island_f1 <- island_f1 +
  theme(legend.position = "none",
        axis.title = element_blank())

site_f1 <- site_f1 +
  theme(legend.position = "none",
        axis.title = element_blank())

# Add a title to each plot using the 'top' argument
island_f1_title <- arrangeGrob(island_f1, top = textGrob("Island scale", 
                                                gp = gpar(fontsize = 13, fontface = "bold")))
site_f1_title <- arrangeGrob(site_f1, top = textGrob("Site scale", 
                                                gp = gpar(fontsize = 13, fontface = "bold")))

grid.arrange(
  arrangeGrob(site_f1_title, island_f1_title, ncol = 2),
  left = textGrob("Balanced accuracy", rot = 90, gp = gpar(fontsize = 13, fontface = "bold"))
)

anova_result <- aov(f1_score ~ layer_comparison, data = result_table)
summary(anova_result)

tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# for ba


## ---- check correlation with similarity in species composition ----


## ---- check correlation with similarity in interaction composition ----

