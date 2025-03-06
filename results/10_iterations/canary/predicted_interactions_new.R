library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

# 1. Create island ID and filter for removed == 1
df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>%
  filter(removed == 1) %>%  mutate(sigm_predicted = sigmoid(predicted_values))

# 2. Aggregate across iterations and island combinations.
#    For each unique interaction (node_from, node_to) compute:
#       - avg_prop: the average proportion of islands where original_links == 1.
#       - avg_sigm_predicted: the average sigm_predicted across all iterations.
# df_summary <- df %>%
#   group_by(node_from, node_to) %>%
#   summarise(
#     avg_prop = mean(original_links == 1, na.rm = TRUE),
#     avg_sigm_predicted = mean(sigm_predicted, na.rm = TRUE),
#     n = n(),  # diagnostic: number of rows contributing to this interaction
#     .groups = "drop"
#   )

# Step 1: For each island and interaction, determine if the interaction was observed.
# Here we use `any(original_links == 1)` so that if the interaction is observed in at least one iteration, we count it.
df_island <- df %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links == 1)),
    # For sigm_predicted, you might take the average across iterations per island.
    island_sigm_predicted = mean(sigm_predicted, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Now, for each unique interaction, compute:
# - The proportion of islands where it was observed.
# - The average predicted probability (averaged over islands).
df_summary <- df_island %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

# 3. Remove interactions that were never observed and predicted as zero.
df_summary <- df_summary %>% 
  filter(!(avg_prop == 0 & avg_sigm_predicted < 0.6)) # needs to be fixed.
# 3. Keep only interactions that are not observed in all islands (avg_prop < 1)
df_summary_not_all <- df_summary %>% filter(avg_prop < 1)

# 4. Order the plants (node_from) by their degree (number of interactions)
#    Here, we count how many unique pollinators each plant has in our aggregated data.
plant_degree <- df_summary_not_all %>%
  group_by(node_from) %>%
  summarise(degree = n(), .groups = "drop") %>%
  arrange(degree)

# Set the order of node_from (plants) as a factor so that the y-axis is ordered by degree.
df_summary_not_all$node_from <- factor(df_summary_not_all$node_from, levels = plant_degree$node_from)

# Optionally, you may also order pollinators (node_to) alphabetically to avoid duplicates.
# df_summary_not_all$node_to <- factor(df_summary_not_all$node_to, levels = sort(unique(df_summary_not_all$node_to)))
pollinator_degree <- df_summary %>%
  group_by(node_to) %>%
  summarise(pollinator_degree = n(), .groups = "drop") %>%
  arrange(desc(pollinator_degree))

# Set node_to factor levels based on pollinator degree (highest degree first)
df_summary$node_to <- factor(df_summary$node_to, levels = pollinator_degree$node_to)
# 5. Plot the heatmap and overlay circles.
#    - The background tiles (geom_tile) show the averaged proportion observed (avg_prop).
#    - We overlay circles (geom_point) for interactions where avg_sigm_predicted > 0.6.
#      Here, we differentiate by color:
#         * Use red if avg_prop == 0 (i.e. never observed empirically)
#         * Use blue if avg_prop > 0
ggplot(df_summary_not_all, aes(x = node_to, y = node_from, fill = avg_prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plant",
    fill = "Prop. observed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
    axis.text.y = element_text(size = 8)
  ) +
  # Overlay circles for high predicted probability (avg_sigm_predicted > 0.6)
  geom_point(
    data = df_summary_not_all %>% filter(avg_sigm_predicted > 0.6, avg_prop == 0),
    aes(x = node_to, y = node_from),
    color = "salmon", shape = 19, size = 3, alpha = 0.8
  ) +
  geom_point(
    data = df_summary_not_all %>% filter(avg_sigm_predicted > 0.6, avg_prop > 0),
    aes(x = node_to, y = node_from),
    color = "thistle", shape = 19, size = 3, alpha = 0.8
  ) +   tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) bquote(italic(.(y[1]) ~ .(y[2]))))) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) bquote(italic(.(y[1]) ~ .(y[2])))))
