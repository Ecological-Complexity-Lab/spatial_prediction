# ---- partner fidelity ----
# here we calculate the "partner fidelity" (measured by Sorensen similarity between the partner composition of each plant and pollinator) and compare it to the accuracy of predictions

# themes
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load libraries ----
library(dplyr)
library(tidyverse)
library(broom)
library(gridExtra)
library(grid)
library(scales)

## ---- load data ----
df <- read_csv('aggregated_equal_0_1_removal_60_1_filtered.csv') # island scale
#df <- read_csv('test_df_itr_1_60_binary.csv') # test on 1 iteration
df <- read_csv('nonbinary_equal_0_1_removal_60_1.csv') # site scale
df <- read_csv('binary_equal_0_1_removal_scaling_island_60_1.csv') # scaled binary version
df <- read_csv('weighted_equal_0_1_removal_scaled_island_60_0.csv') # scaled weighted version


df <- df %>% filter(k == 2 & lambda == 0.1) # for site scale

# filter out cases in which train = test layer
df_fidelity <- df %>% filter(train_layer == test_layer) %>% 
  filter(original_links == 1) %>% filter(itr == 1)

## ---- plants fidelity ----
# 1. For each plant (node_from) and layer, gather the pollinators (node_to).
#    (Assuming 'train_layer' is the relevant layer ID—adapt as needed.)
df_plant_partners <- df_fidelity %>%
  distinct(node_from, train_layer, node_to) %>%
  group_by(node_from, train_layer) %>%
  summarise(partners = list(unique(node_to)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen <- df_plant_partners %>%
  group_by(node_from) %>%
  summarise(
    mean_sorensen_plants = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

df_merged <- df %>%
  left_join(df_sorensen, by = "node_from")

## ---- pollinators fidelity ----

df_pollinators <- df_fidelity %>%
  distinct(node_to, train_layer, node_from) %>%
  group_by(node_to, train_layer) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen_pollinators <- df_pollinators %>%
  group_by(node_to) %>%
  summarise(
    mean_sorensen_pollinators = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

df_merged <- df_merged %>%
  left_join(df_sorensen_pollinators, by = "node_to")

write.csv(df_sorensen_pollinators, "sorensen_pollinators_island.csv")
write.csv(df_sorensen, "sorensen_plants_island.csv")

# weighted version
## ---- add distance between predicted and observed values ----
df_merged <- df_merged %>% mutate(predicted_value_sigm = sigmoid(predicted_values))
df_merged <- df_merged %>% mutate(predicted_sigm_obs = predicted_value_sigm - original_links)
df_merged <- df_merged %>% mutate(pred_obs_abs = abs(predicted_sigm_obs))

df_merged <- df_merged %>% mutate(pred_obs_abs = abs(predicted_values - original_links))
df_merged_removed <- df_merged %>% filter(removed == 1)

# here you can separate diagonals and off-diagonals
df_merged_removed_offs <- df_merged_removed %>% filter(train_layer != test_layer)
df_merged_removed_offs_itr1 <- df_merged_removed_offs %>% filter(itr == 1)
df_merged_removed_diag <- df_merged_removed %>% filter(train_layer == test_layer)

## ---- correlation pollinators ----
#nonbinary_fidelity_merged_removed_itr_1 <- nonbinary_fidelity_merged %>% filter(itr == 1) # try 1 itr
correlation <- cor.test(df_merged_removed$pred_obs_abs, df_merged_removed$mean_sorensen_pollinators, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

pollinator_fidelity_cor <- ggplot(df_merged_removed, aes(x = mean_sorensen_pollinators, y = pred_obs_abs)) +
  geom_point(color = "thistle", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "Predicted - observed") +
       #title = "Island scale (pollinators) - diag all itr") +
  tme +
  annotate("text",
           x = 0.77, y = 75,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black")

## ---- correlation plants ----
#nonbinary_fidelity_merged_removed_itr_1 <- nonbinary_fidelity_merged %>% filter(itr == 1) # try 1 itr
correlation <- cor.test(df_merged_removed$pred_obs_abs, df_merged_removed$mean_sorensen_plants, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

plant_fidelity_cor <- ggplot(df_merged_removed, aes(x = mean_sorensen_plants, y = pred_obs_abs)) +
  geom_point(color = "darkseagreen3", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "Predicted - observed") +
       #title = "Island scale (plants) - diag all itr") +
  tme +
  annotate("text",
           x = 0.12, y = 75,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black") +
  theme(axis.title.y = element_blank()) # for the unified plot


# Combine the plots and add a shared y-axis title
combined_plot <- plant_fidelity_cor + pollinator_fidelity_cor + 
  plot_layout(guides = "collect") & 
  theme(axis.title.y = element_text(size = 14))  # Adjust y-axis title style

combined_plot
# Add a common y-axis title
combined_plot + plot_annotation(title = "Shared Y-axis Title")

## ---- check if the correlation is the same for the different iterations ----

# Assuming your dataframe is called df
cor_results <- df_merged_removed %>%
  group_by(train_layer, test_layer, itr) %>%
  summarise(
    cor_test = list(cor.test(pred_obs_abs, mean_sorensen_pollinators, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(
    correlation = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  select(train_layer, test_layer, itr, correlation, p_value)

# Print the results
print(cor_results)
write.csv(cor_results, 'pollinator_fidelity_itr_corr.csv')

## ---- correlate fidelity with evaluators ----

summary_df <- df_merged %>%
  group_by(train_layer, test_layer) %>%
  summarise(
    avg_sorensen_plants = mean(mean_sorensen_plants, na.rm = TRUE),
    avg_sorensen_pollinators = mean(mean_sorensen_pollinators, na.rm = TRUE)
  ) %>%
  ungroup()

# Join the summarized data with working_df based on train_layer and test_layer
working_df <- read.csv("working_df_islands_evaluators_distance.csv")
working_df <- read.csv("result_canaries_distance_names_jaccard_island_scaled.csv") # scaled binary version
working_df <- read.csv("working_df_island_weighted_scaled_evaluators_distance.csv") # scaled weighted version

working_df <- working_df %>%
  left_join(summary_df, by = c("train_layer", "test_layer"))

# View the updated working_df
view(working_df)
write.csv(working_df, "working_df_weighted_scaled_island_distance_fidelity.csv")

working_df_offs <- working_df %>% filter (train_layer != test_layer)
correlation <- cor.test(working_df_offs$f1_score, working_df_offs$avg_sorensen_pollinators, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

pollinator_fidelity_cor <- ggplot(working_df_offs, aes(x = avg_sorensen_pollinators, y = f1_score)) +
  geom_point(color = "thistle", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "F1 score",
       title = "Island scale (pollinators)") +
  tme +
  annotate("text",
           x = 0.31, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black")+
  theme(axis.title.y = element_blank()) # for the unified plot

correlation <- cor.test(working_df_offs$f1_score, working_df_offs$avg_sorensen_plants, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

plant_fidelity_cor <- ggplot(working_df_offs, aes(x = avg_sorensen_plants, y = f1_score)) +
  geom_point(color = "darkseagreen3", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "F1 score",
       title = "Island scale (plants)") +
  tme +
  annotate("text",
           x = 0.17, y = 0.7,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black") +
  theme(axis.title.y = element_blank()) # for the unified plot


grid.arrange(
  arrangeGrob(plant_fidelity_cor, pollinator_fidelity_cor, ncol = 2),
  left = textGrob("F1 score", rot = 90, gp = gpar(fontsize = 13, fontface = "bold"))
)
