# ---- comparing scales ----
## ---- load libraries ----
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

## ---- functions ----
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
## ---- load data ----
site_scale <- read_csv('result_summary_canary_with_distance.csv') # site scale
d <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv') # prepare evaluators for island scale

d <- d %>%
  filter(removed == 1) %>% 
  filter(k == 2) %>% 
  filter(lambda == 0.1) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0))

island_scale <- d %>%
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
  ungroup() %>% 
  write_csv('evaluation_df_all_itr_60_binary_island.csv')

# add distances


## ---- summarize predictive performance ----
# summarize max, min, average, se of f1, ba, recall for both scales
island_scale_off <- island_scale %>% filter(test_layer != train_layer) # these are for off-diagonals. for all data use site_scale and island_scale
site_scale_off <- site_scale %>% filter(test_layer != train_layer)

# for diagonal:
island_scale_diag <- island_scale %>% filter(test_layer == train_layer) # these are for diagonals. for all data use site_scale and island_scale
site_scale_diag <- site_scale %>% filter(test_layer == train_layer)

island_scale <- read.csv('evaluation_df_all_itr_60_binary_island.csv')

df1_labeled <- site_scale %>%
  mutate(scale = "Site")

df2_labeled <- island_scale %>%
  mutate(scale = "Island")

df_combined <- bind_rows(df1_labeled, df2_labeled)

df_long <- df_combined %>%
  pivot_longer(
    cols = c("f1_score", "recall", "precision", "balanced_accuracy", "mcc"),
    names_to = "metric",
    values_to = "value"
  )


ggplot(df_long, aes(x = metric, y = value, fill = scale)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    title = "Comparison of Performance (Off-diags + Diag)",
    x = "Metric",
    y = "Value"
  ) +
  scale_fill_manual(values = c("Site" = "lightsteelblue2", "Island" = "wheat2")) +  # Custom colors
  scale_x_discrete(labels = c(
    "f1_score" = "F1 fcore",
    "recall" = "Recall",
    "precision" = "Precision",
    "mcc" = "MCC",
    "balanced_accuracy" = "Balanced \naccuracy"
  )) +
  tme

## ---- test for significance ----

metrics <- c("f1_score", "recall", "precision", "balanced_accuracy", "mcc")

results <- lapply(metrics, function(metric) {
  test_normality_site <- shapiro.test(site_scale[[metric]])$p.value
  test_normality_island <- shapiro.test(island_scale[[metric]])$p.value
  
  if (test_normality_site > 0.05 & test_normality_island > 0.05) {
    test <- t.test(site_scale[[metric]], island_scale[[metric]], var.equal = FALSE)
  } else {
    test <- wilcox.test(site_scale[[metric]], island_scale[[metric]])
  }
  
  data.frame(
    Metric = metric,
    Test = ifelse(test_normality_site > 0.05 & test_normality_island > 0.05, "T-test", "Wilcoxon"),
    P_value = test$p.value
  )
})

results_df <- do.call(rbind, results)
print(results_df)

# Prepare dataset with labels

## plot with significance levels

# Prepare dataset with labels
df1_labeled <- site_scale %>% mutate(scale = "Site")
df2_labeled <- island_scale %>% mutate(scale = "Island")
df_combined <- bind_rows(df1_labeled, df2_labeled)

df_long <- df_combined %>%
  pivot_longer(cols = c("f1_score", "recall", "precision", "balanced_accuracy", "mcc"),
               names_to = "metric",
               values_to = "value")

# Define significance function
get_pvalue_asterisks <- function(p) {
  if (p < 0.001) return("***")  # Highly significant
  else if (p < 0.01) return("**")  # Very significant
  else if (p < 0.05) return("*")  # Significant
  else return("ns")  # Not significant
}

# Perform statistical tests and collect results
metrics <- c("f1_score", "recall", "precision", "balanced_accuracy", "mcc")

stat_results <- lapply(metrics, function(metric) {
  data_metric <- df_long %>% filter(metric == !!metric)  # Filter for the specific metric
  
  test <- t.test(value ~ scale, data = data_metric)  # Perform t-test
  
  p_value <- test$p.value
  significance <- get_pvalue_asterisks(p_value)
  
  data.frame(
    metric = metric,
    p_value = p_value,
    significance = significance
  )
})

stat_results_df <- do.call(rbind, stat_results)

# Merge significance levels with the dataset
df_long <- df_long %>%
  left_join(stat_results_df, by = "metric")

# Create the boxplot with significance annotations
ggplot(df_long, aes(x = metric, y = value, fill = scale)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Comparison of Performance", x = "Metric", y = "Value") +
  scale_fill_manual(values = c("Site" = "lightsteelblue2", "Island" = "wheat2")) +  # Custom colors
  scale_x_discrete(labels = c(
    "f1_score" = "F1 score",
    "recall" = "Recall",
    "precision" = "Precision",
    "balanced_accuracy" = "Balanced \naccuracy",
    "mcc" = "MCC"
  )) +  # Properly formatted labels
  stat_compare_means(aes(group = scale), method = "t.test", label = "p.signif", 
                     label.y = max(df_long$value, na.rm = TRUE) + 0.05,
                     size = 5)  + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14), ) + tme

library(ggplot2)
library(ggpubr)  # For stat_compare_means()

ggplot(df_long, aes(x = metric, y = value, fill = scale)) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "Comparison of Performance", x = "Metric", y = "Value") +
  scale_fill_manual(values = c("Site" = "lightsteelblue2", "Island" = "wheat2")) +  # Custom colors
  scale_x_discrete(labels = c(
    "f1_score" = "F1 score",
    "recall" = "Recall",
    "precision" = "Precision",
    "balanced_accuracy" = "Balanced \naccuracy",
    "mcc" = "MCC"
  )) +  # Properly formatted labels
  stat_compare_means(
    aes(group = scale), method = "t.test", label = "p.signif", 
    label.y = max(df_long$value, na.rm = TRUE) + 0.05, 
    size = 8  # Increase significance asterisk size
  ) +
  theme(
    legend.text = element_text(size = 14)  # Increase legend font size
  )

## ---- heatmap of island scale ----
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Display the updated tibble
print(net_name)

# Extract island/subregion names (remove "site X")
net_name_island <- net_name %>%
  mutate(name = gsub(" site [12]", "", name)) %>%  # Remove " site 1" and " site 2"
  distinct(name) %>%  # Keep unique names
  mutate(layer_id = row_number())  # Assign new layer_id

# Print the cleaned list
print(net_name_island)

island_scale <- island_scale %>%
  mutate(train_layer = as.integer(as.character(train_layer)),  # Ensure train_layer is an integer
         test_layer = as.integer(as.character(test_layer))) %>% # Ensure test_layer is an integer
  left_join(net_name_island, by = c("train_layer" = "layer_id")) %>% 
  rename(train_layer_name = name) %>%                          # Use a different name for clarity
  left_join(net_name_island, by = c("test_layer" = "layer_id")) %>%  
  rename(test_layer_name = name)    

#write.csv(island_scale, "working_df_evaluators_island_names.csv")
island_scale$diagonal <- island_scale$train_layer == island_scale$test_layer

layer_to_layer_plot_islands <- 
  ggplot(island_scale, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = island_scale[island_scale$train_layer == island_scale$test_layer, ],
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

print(layer_to_layer_plot_islands)

## ---- correlate performance with distance ----
# first convert distances to distances between islands
island_scale <- read.csv('working_df_evaluators_island_names.csv')
distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL)

# Function to extract island names (removes "_site_X")
extract_island <- function(name) {
  gsub("_site_[12]", "", name)
}

# Create new table with averaged distances at the island level
distance_island_table <- distance_table %>%
  mutate(
    from_island = extract_island(from),
    to_island = extract_island(to)
  ) %>%
  group_by(from_island, to_island) %>%
  summarise(
    avg_distance_m = mean(distance_m),
    avg_distance_km = mean(distance_km),
    .groups = "drop"
  ) %>%
  mutate(
    avg_distance_m = ifelse(from_island == to_island, 0, avg_distance_m),
    avg_distance_km = ifelse(from_island == to_island, 0, avg_distance_km)
  ) %>%
  rename(from = from_island, to = to_island)  # Rename after calculation

# Print result
print(distance_island_table)

# Modify the 'from' and 'to' columns in distance_island_table
distance_island_table <- distance_island_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

# Perform the left join
island_scale <- island_scale %>%
  left_join(distance_island_table, by = c("train_layer_name" = "from", "test_layer_name" = "to")) %>%
  mutate(distance_km = avg_distance_km)

island_scale <- island_scale %>% select(-avg_distance_km)

write.csv(island_scale, 'working_df_islands_evaluators_distance.csv')

correlation <- cor.test(island_scale$f1_score, island_scale$distance_km, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, format = "f", digits = 3)

#p_value <- round(correlation$p.value, 5)

label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot_canary <- 
  ggplot(island_scale, aes(x = distance_km, y = f1_score)) +
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

island_scale_off <- island_scale %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer != test_layer)

correlation <- cor.test(island_scale_off$f1_score, island_scale_off$distance_km, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, format = "f", digits = 3)

#p_value <- round(correlation$p.value, 5)

label_text <- paste0("r = ", r_value, ", p = ", p_value)

# Plot the correlation between balanced accuracy and geographic distance

cor_plot_canary_off <- 
  ggplot(island_scale_off, aes(x = distance_km, y = f1_score)) +
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

print(cor_plot_canary_off)

## ---- distribution of evaluators ----
site_df <- read.csv('result_summary_canary_with_distance.csv')
result_summary_isl <- read.csv('working_df_island_distance_fidelity_jaccard_netsize.csv')

site_distrib <- ggplot(site_df, aes(x = balanced_accuracy)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.5) + 
  labs(title = "Site scale",
       x = "Balanced Accuracy",
       y = "Count") +
  tme

isl_distrib <- ggplot(result_summary_isl, aes(x = recall)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.5) + 
  labs(title = "Island scale",
       x = "Recall",
       y = "Count") +
  tme

