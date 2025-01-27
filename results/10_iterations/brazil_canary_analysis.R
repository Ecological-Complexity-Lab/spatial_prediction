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
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) 

d %>% 
  ggplot(aes(original_links,predicted_prob_sigm))+geom_point()

## ---- some evaluators ----


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
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) 

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
# Rename joined column

# Set factor levels for training and predicting layers
layer_levels <- as.character(1:7) # for Brazil
layer_levels <- as.character(1:14) # for Canary Islands
result_summary$train_layer <- factor(result_summary$train_layer, levels = layer_levels)
result_summary$test_layer <- factor(result_summary$test_layer, levels = unique(result_summary$test_layer))


result_summary$diagonal <- result_summary$train_layer == result_summary$test_layer

layer_to_layer_plot_all_itr_brzail <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
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

print(layer_to_layer_plot_all_itr_brzail)

layer_to_layer_plot_all_itr_canary <- 
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

print(layer_to_layer_plot_all_itr_canary)


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
result_summary[1:15, c("train_layer","test_layer","geographic_distance")]

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
  labs(x = "Geographic distance (km)",
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
  labs(x = "Geographic distance (km)",
       y = "F1 score") +
  tme + 
  annotate("text",
           x = 7, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")

print(cor_plot)

# trying the canary islands...

# Load the pairwise distance matrix
distance_matrix <- read.csv("distance_between_sites_canary.csv", row.names = NULL) # we need to match the layer names to the numbering in the results table
# # acast syntax: acast(data, row_variable ~ col_variable, value.var = "...")
distance_matrix <- acast(distance_matrix, from ~ to, value.var = "distance_km")
distance_matrix[is.na(dist_mat)] <- 0
diag(distance_matrix) <- 0

# For each row's train_layer ID, find which row in net_name has the same layer_id,
# and pull out the corresponding 'name'.
result_summary$train_name <- net_name$name[ match(result_summary$train_layer, net_name$layer_id) ]

# Same for test_layer
result_summary$test_name <- net_name$name[ match(result_summary$test_layer, net_name$layer_id) ]

library(dplyr)

result_summary <- result_summary %>%
  rowwise() %>%                                # loop over rows
  mutate(
    distance_km = distance_matrix[train_name,  # train_name is e.g. "El_Hierro_site_1"
                                  test_name]   # test_name is e.g. "Gran_Canaria_site_2"
  ) %>%
  ungroup()


# # site_map$layer_id are the integer IDs
# # site_map$name are the site names
# site_map_vector <- setNames(net_name$layer_id, net_name$name)
# # Rename rows
# rownames(distance_matrix) <- site_map_vector[ rownames(distance_matrix) ]
# # Rename columns
# colnames(distance_matrix) <- site_map_vector[ colnames(distance_matrix) ]
# diag(distance_matrix) <- 0


distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- colnames(distance_matrix) <- 1:14

result_summary$train_layer <- as.factor(result_summary$train_layer)
result_summary$test_layer <- as.factor(result_summary$test_layer)

# Create a column with the geographic distance between the layers for each row
result_summary <- result_summary %>%
  mutate(geographic_distance = mapply(function(train, test) distance_matrix[as.character(train), as.character(test)],
                                      train_layer, test_layer))
result_summary[1:15, c("train_layer","test_layer","geographic_distance")]

correlation <- cor.test(result_summary$balanced_accuracy, result_summary$geographic_distance, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot <- 
  ggplot(result_summary, aes(x = geographic_distance, y = balanced_accuracy)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographic distance (km)",
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
  labs(x = "Geographic distance (km)",
       y = "F1 score") +
  tme + 
  annotate("text",
           x = 7, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")

print(cor_plot)
