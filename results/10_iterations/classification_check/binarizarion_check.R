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

# Min-max normalization function
normalize_min_max <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

## ---- checking binarization on 1 iteration ----

# d <- read_csv('nonbinary_equal_0_1_removal_25_1.csv') # brazil
# d <- read_csv('nonbinary_equal_0_1_removal_60_1.csv') # canary islands
d <- read_csv('test_df_itr_1_60_binary.csv') # canary islands test df, contains 1 iteration, k==2, lambda == 0.01

d <- d %>%
  filter(removed == 1) %>% 
  # filter(k == 2) %>% 
  # filter(lambda == 0.1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) %>% 
  mutate(predicted_prob_min_max = normalize_min_max(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_minmax = if_else(predicted_values > 0.5, 1, 0)) %>% 
  mutate(predicted_prob_above_0 = if_else(predicted_values > 0, 1, 0))  %>% # positive predicted values are ones
  mutate(predicted_bin_sigm0.5 = ifelse(predicted_prob_sigm == 0.5, 0, 1)) # similarly to the classification in the older version code, sigmoid values of 0.5 are zeros

result_summary <- d %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    # Metrics for predicted_bin_sigm
    TP_sigm = sum(original_links == 1 & predicted_bin_sigm == 1),
    FN_sigm = sum(original_links == 1 & predicted_bin_sigm == 0),
    TN_sigm = sum(original_links == 0 & predicted_bin_sigm == 0),
    FP_sigm = sum(original_links == 0 & predicted_bin_sigm == 1),
    specificity_sigm = TN_sigm / (TN_sigm + FP_sigm),
    precision_sigm = TP_sigm / (TP_sigm + FP_sigm),
    recall_sigm = TP_sigm / (TP_sigm + FN_sigm),
    f1_score_sigm = 2 * (precision_sigm * recall_sigm) / (precision_sigm + recall_sigm),
    balanced_accuracy_sigm = (recall_sigm + specificity_sigm) / 2,
    mcc_sigm = (TP_sigm * TN_sigm - FP_sigm * FN_sigm) / sqrt((TP_sigm + FP_sigm) * (TP_sigm + FN_sigm) * (TN_sigm + FP_sigm) * (TN_sigm + FN_sigm)),
    
    # Metrics for predicted_bin_minmax
    TP_minmax = sum(original_links == 1 & predicted_bin_minmax == 1),
    FN_minmax = sum(original_links == 1 & predicted_bin_minmax == 0),
    TN_minmax = sum(original_links == 0 & predicted_bin_minmax == 0),
    FP_minmax = sum(original_links == 0 & predicted_bin_minmax == 1),
    specificity_minmax = TN_minmax / (TN_minmax + FP_minmax),
    precision_minmax = TP_minmax / (TP_minmax + FP_minmax),
    recall_minmax = TP_minmax / (TP_minmax + FN_minmax),
    f1_score_minmax = 2 * (precision_minmax * recall_minmax) / (precision_minmax + recall_minmax),
    balanced_accuracy_minmax = (recall_minmax + specificity_minmax) / 2,
    mcc_minmax = (TP_minmax * TN_minmax - FP_minmax * FN_minmax) / sqrt((TP_minmax + FP_minmax) * (TP_minmax + FN_minmax) * (TN_minmax + FP_minmax) * (TN_minmax + FN_minmax)),
    
    # Metrics for predicted_prob_above_0
    TP_prob0 = sum(original_links == 1 & predicted_prob_above_0 == 1),
    FN_prob0 = sum(original_links == 1 & predicted_prob_above_0 == 0),
    TN_prob0 = sum(original_links == 0 & predicted_prob_above_0 == 0),
    FP_prob0 = sum(original_links == 0 & predicted_prob_above_0 == 1),
    specificity_prob0 = TN_prob0 / (TN_prob0 + FP_prob0),
    precision_prob0 = TP_prob0 / (TP_prob0 + FP_prob0),
    recall_prob0 = TP_prob0 / (TP_prob0 + FN_prob0),
    f1_score_prob0 = 2 * (precision_prob0 * recall_prob0) / (precision_prob0 + recall_prob0),
    balanced_accuracy_prob0 = (recall_prob0 + specificity_prob0) / 2,
    mcc_prob0 = (TP_prob0 * TN_prob0 - FP_prob0 * FN_prob0) / sqrt((TP_prob0 + FP_prob0) * (TP_prob0 + FN_prob0) * (TN_prob0 + FP_prob0) * (TN_prob0 + FN_prob0)),
    
    # Metrics for predicted_bin_sigm0.5
    TP_sigm05 = sum(original_links == 1 & predicted_bin_sigm0.5 == 1),
    FN_sigm05 = sum(original_links == 1 & predicted_bin_sigm0.5 == 0),
    TN_sigm05 = sum(original_links == 0 & predicted_bin_sigm0.5 == 0),
    FP_sigm05 = sum(original_links == 0 & predicted_bin_sigm0.5 == 1),
    specificity_sigm05 = TN_sigm05 / (TN_sigm05 + FP_sigm05),
    precision_sigm05 = TP_sigm05 / (TP_sigm05 + FP_sigm05),
    recall_sigm05 = TP_sigm05 / (TP_sigm05 + FN_sigm05),
    f1_score_sigm05 = 2 * (precision_sigm05 * recall_sigm05) / (precision_sigm05 + recall_sigm05),
    balanced_accuracy_sigm05 = (recall_sigm05 + specificity_sigm05) / 2,
    mcc_sigm05 = (TP_sigm05 * TN_sigm05 - FP_sigm05 * FN_sigm05) / sqrt((TP_sigm05 + FP_sigm05) * (TP_sigm05 + FN_sigm05) * (TN_sigm05 + FP_sigm05) * (TN_sigm05 + FN_sigm05))
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    across(starts_with("specificity"), mean, na.rm = TRUE),
    across(starts_with("precision"), mean, na.rm = TRUE),
    across(starts_with("recall"), mean, na.rm = TRUE),
    across(starts_with("f1_score"), mean, na.rm = TRUE),
    across(starts_with("balanced_accuracy"), mean, na.rm = TRUE),
    across(starts_with("mcc"), mean, na.rm = TRUE)
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

# # Set factor levels for training and predicting layers
# layer_levels <- as.character(1:7) # for Brazil
# layer_levels <- as.character(1:14) # for Canary Islands
# result_summary$train_layer <- factor(result_summary$train_layer, levels = layer_levels)
# result_summary$test_layer <- factor(result_summary$test_layer, levels = unique(result_summary$test_layer))

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

layer_to_layer_plot_canary_sigm <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f1_score_sigm)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
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

print(layer_to_layer_plot_canary_sigm)

layer_to_layer_plot_canary_minmax <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f1_score_minmax)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
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

print(layer_to_layer_plot_canary_minmax)

layer_to_layer_plot_canary_prob0 <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f1_score_prob0)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
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

layer_to_layer_plot_canary_sigm05 <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = f1_score_sigm05)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
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

print(layer_to_layer_plot_all_itr_canary)

library(gridExtra)

grid.arrange(
  layer_to_layer_plot_canary_sigm + ggtitle("Sigmoid"),
  layer_to_layer_plot_canary_minmax + ggtitle("Min-Max"),
  layer_to_layer_plot_canary_prob0 + ggtitle("Prob > 0"),
  layer_to_layer_plot_canary_sigm05 + ggtitle("Sigmoid 0.5"),
  ncol = 2,
  nrow = 2 # Arrange in a single row
)


## ---- check if the differences are significant ----

result_summary <- result_summary %>%
  mutate(layer_comparison = ifelse(train_layer == test_layer, "Diagonal", "Off-Diagonal"))

library(ggpubr)

# Define a function to create boxplots for each metric
plot_f1_boxplot <- function(metric) {
  ggplot(result_summary, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Comparison of", metric, "by Layer Type"),
         x = "Layer Comparison",
         y = metric) +
    theme(legend.position = "none") +  # Hide redundant legend
    stat_compare_means(method = "t.test", label = "p.signif") +  # Adds p-value annotation
    tme
    }

# Plot for f1_score_sigm
boxplot_sigm <- plot_f1_boxplot("f1_score_sigm")

# Plot for f1_score_prob0
boxplot_prob0 <- plot_f1_boxplot("f1_score_prob0")

# Plot for f1_score_sigm05
boxplot_sigm05 <- plot_f1_boxplot("f1_score_sigm05")

grid.arrange(
  boxplot_sigm + ggtitle("Sigmoid"),
  boxplot_prob0 + ggtitle("Prob > 0"),
  boxplot_sigm05 + ggtitle("Sigmoid 0.5"),
  ncol = 3
)

# Function to perform statistical tests
compare_f1_scores <- function(metric) {
  t_test <- t.test(result_table[[metric]] ~ result_table$layer_comparison)
  wilcox_test <- wilcox.test(result_table[[metric]] ~ result_table$layer_comparison)
  
  cat("Results for", metric, "\n")
  cat("T-test p-value:", t_test$p.value, "\n")
  cat("Wilcoxon test p-value:", wilcox_test$p.value, "\n\n")
}

# Compare F1 scores statistically
compare_f1_scores("f1_score_sigm")
compare_f1_scores("f1_score_prob0")
compare_f1_scores("f1_score_sigm05")


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
result_summary %>% 
  write_csv('result_summary_canary_with_distance.csv') # summary of all iterations

# Load the pairwise distance matrix
distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL) # we need to match the layer names to the numbering in the results table

# For each row's train_layer ID, find which row in net_name has the same layer_id,
# and pull out the corresponding 'name'.
result_summary$train_name <- net_name$name[ match(result_summary$train_layer, net_name$layer_id) ]

# Same for test_layer
result_summary$test_name <- net_name$name[ match(result_summary$test_layer, net_name$layer_id) ]

distance_table_sym <- distance_table %>%
  bind_rows(
    distance_table %>%
      rename(to = from, from = to)  # swap columns
  ) %>%
  distinct(from, to, distance_km)

# result_summary has train_name, test_name, 
# plus an old "distance_km" you want to replace

result_summary <- result_summary %>%
  select(-distance_km) %>%                # remove old distance_km if it exists
  left_join(
    distance_table_sym,
    by = c("train_name" = "from", "test_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_name == test_name,
                               0,              # distance = 0 if same site
                               distance_km))   # otherwise, keep joined distance



result_summary %>%
  anti_join(distance_table_sym, by = c("train_name" = "from", "test_name" = "to"))

correlation <- cor.test(result_summary$f1_score, result_summary$distance_km, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 5)  # or round as you prefer
p_value <- formatC(correlation$p.value, format = "f", digits = 5)

p_value <- round(correlation$p.value, 3)

label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot_canary <- 
  ggplot(result_summary, aes(x = distance_km, y = f1_score)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographic distance (km)",
       y = "F1 score") +
  annotate("text",
           x = 370, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black") +
  tme

print(cor_plot_canary)

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
