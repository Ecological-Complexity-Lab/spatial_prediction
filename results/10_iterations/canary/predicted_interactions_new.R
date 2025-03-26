library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(gridExtra)

source("~/Documents/github/softimpute/results/useful_for_plotting.R")

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 12, color = "black"),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))

df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

# 1. Create island ID and filter for removed == 1
df <- df %>%
  filter(removed == 1) %>% mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(island_id = paste(train_layer, test_layer, sep = "_"))

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
  filter(!(avg_prop == 0 & avg_sigm_predicted < 0.5))
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
    fill = "Proportion\nof islands\nobserved"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 8) 
  ) + 
  # Overlay circles for high predicted probability (avg_sigm_predicted > 0.6)
  geom_point(
    data = df_summary_not_all %>% filter(avg_sigm_predicted > 0.5, avg_prop == 0),
    aes(x = node_to, y = node_from),
    color = "thistle", shape = 19, size = 3, alpha = 0.6
  ) +
  geom_point(
    data = df_summary_not_all %>% filter(avg_sigm_predicted > 0.5, avg_prop > 0),
    aes(x = node_to, y = node_from),
    color = "salmon", shape = 1, size = 4.5, stroke = 1, alpha = 0.7
  ) +   tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))
 
## ---- always true interactions and relation to degree ----
### ---- load data ----
df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

# 1. Create island ID and filter for removed == 1
df <- df %>%
  mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(island_id = paste(train_layer, test_layer, sep = "_"))

### ---- calculate degree ----
# Keep only the rows where train_layer == test_layer (if that’s what defines an island)
df_isalnd <- df %>% filter(train_layer == test_layer)

# For plants: Calculate island-level degree (number of pollinators with original_links == 1)
plant_degree_island <- df_isalnd %>%
  filter(original_links == 1) %>%
  group_by(island_id, node_from) %>%
  summarise(degree = n_distinct(node_to), .groups = "drop")

# Then average the degree per plant across islands
plant_avg_degree <- plant_degree_island %>%
  group_by(node_from) %>%
  summarise(avg_degree = mean(degree), .groups = "drop")

# For pollinators: Calculate island-level degree (number of plants with original_links == 1)
poll_degree_island <- df_isalnd %>%
  filter(original_links == 1) %>%
  group_by(island_id, node_to) %>%
  summarise(degree = n_distinct(node_from), .groups = "drop")

# Then average the degree per pollinator across islands
poll_avg_degree <- poll_degree_island %>%
  group_by(node_to) %>%
  summarise(avg_degree = mean(degree), .groups = "drop")

plants_long <- plant_avg_degree %>%
  rename(species = node_from) %>%
  mutate(role = "plant")

polls_long <- poll_avg_degree %>%
  rename(species = node_to) %>%
  mutate(role = "pollinator")

df_species_avg_degree <- bind_rows(plants_long, polls_long)

### ---- identify always correct predictions ----
# First, check prediction correctness per iteration within each island
df_iter <- df %>%
  group_by(node_from, node_to, island_id, itr) %>%
  summarise(
    iter_correct = as.integer((original_links == 1 & sigm_predicted > 0.5) |
                                (original_links == 0 & sigm_predicted <= 0.5)),
    .groups = "drop"
  )

# For each island, the interaction is correct if every iteration was correct
df_island_t <- df_iter %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    island_correct = as.integer(all(iter_correct == 1)),
    .groups = "drop"
  )

# Finally, an interaction is always correct if it is correct in every island
df_interaction_summary <- df_island_t %>%
  group_by(node_from, node_to) %>%
  summarise(
    always_correct = as.integer(all(island_correct == 1)),
    prop_correct = mean(island_correct),  # proportion of islands correct (for information)
    n_islands = n(),
    .groups = "drop"
  )

# For plants:
df_plants_correct <- df_interaction_summary %>%
  group_by(node_from) %>%
  summarise(
    prop_always_correct = mean(always_correct),
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  rename(species = node_from) %>%
  mutate(role = "plant")

# For pollinators:
df_pollinators_correct <- df_interaction_summary %>%
  group_by(node_to) %>%
  summarise(
    prop_always_correct = mean(always_correct),
    n_interactions = n(),
    .groups = "drop"
  ) %>%
  rename(species = node_to) %>%
  mutate(role = "pollinator")

df_correct_species <- bind_rows(df_plants_correct, df_pollinators_correct)

df_species_final <- df_correct_species %>%
  left_join(df_species_avg_degree, by = c("species", "role"))

### ---- plot ----
ggplot(df_species_final, aes(x = avg_degree, y = prop_always_correct, color = role)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ role) +
  labs(
    x = "Average Degree per Island",
    y = "Proportion of Interactions Always Correct",
    title = "Relationship Between Island-Averaged Degree and Prediction Accuracy",
    color = "Role"
  ) +
  theme_minimal()

## ---- obs - pred = 0 ----

# Step 1: Create the new columns based on the conditions
df <- df %>%
  mutate(
    TP = ifelse(original_links == 1 & sigm_predicted > 0.5, 1, 0),
    TN = ifelse(original_links == 0 & sigm_predicted <= 0.5, 1, 0),
    FP = ifelse(original_links == 0 & sigm_predicted > 0.5, 1, 0),
    FN = ifelse(original_links == 1 & sigm_predicted <= 0.5, 1, 0)
  )

# Step 2: Group by the relevant columns and summarize across iterations
layer_summary <- df %>%
  group_by(train_layer, test_layer, node_to, node_from) %>%
  summarise(
    count_itr = n(),
    sum_TP = sum(TP),
    sum_TN = sum(TN),
    sum_FP = sum(FP),
    sum_FN = sum(FN)
  ) %>%
  ungroup() %>%
  mutate(
    consistent_TP = (count_itr == sum_TP),
    consistent_TN = (count_itr == sum_TN)
  )

# which interactions are always predicted accurately?
final_summary <- layer_summary %>%
  group_by(node_to, node_from) %>%
  summarise(
    all_consistent_TP = all(consistent_TP),
    all_consistent_TN = all(consistent_TN)
  )

print(final_summary)

# View the summarized data frame
summary_df_TP <- final_summary %>% filter(all_consistent_TP == TRUE)
summary_df_TN <- final_summary %>% filter(all_consistent_TN == TRUE)

# Step 1: (Assuming you already have computed consistency per layer combination)
# For example, layer_summary has columns: train_layer, test_layer, node_to, node_from, consistent_TP, consistent_TN.
# Here, we flag an interaction as "accurate" if it's either consistently TP or TN:
layer_summary <- layer_summary %>%
  mutate(accurate_interaction = consistent_TP | consistent_TN)

# Step 2: Create a table of species that appear in any accurate interaction per layer combination
species_accurate <- layer_summary %>%
  filter(accurate_interaction) %>%
  select(train_layer, test_layer, node_to, node_from) %>%
  pivot_longer(cols = c(node_to, node_from), 
               names_to = "role", 
               values_to = "species") %>%
  distinct(train_layer, test_layer, species) %>%
  mutate(accurate = TRUE)

# Step 3: For degree calculation, get the unique network (collapse 10 iterations per layer combination)
distinct_network <- df %>%
  group_by(train_layer, test_layer, node_to, node_from) %>%
  summarise(original_links = first(original_links), .groups = "drop") %>%
  filter(original_links == 1)

# Compute degree for each species within each layer combination.
# Calculate separately for when species is in node_to and node_from, then combine.
degree_to <- distinct_network %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(degree = n(), .groups = "drop") %>%
  rename(species = node_to)

degree_from <- distinct_network %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(degree = n(), .groups = "drop") %>%
  rename(species = node_from)

degree_all <- bind_rows(degree_to, degree_from) %>%
  group_by(train_layer, test_layer, species) %>%
  summarise(degree = sum(degree), .groups = "drop")

# Step 4: Merge degree information with the accurate involvement flag
analysis_df <- degree_all %>%
  left_join(species_accurate, by = c("train_layer", "test_layer", "species")) %>%
  # species not flagged get NA, so replace with FALSE
  mutate(accurate = ifelse(is.na(accurate), FALSE, accurate))

# Step 5: Compare the degree distributions
# For example, using a boxplot:
ggplot(analysis_df, aes(x = as.factor(accurate), y = degree)) +
  geom_point() +
  labs(x = "Species Involved in Accurate Interactions (TP or TN)",
       y = "Degree (Number of Interactions)",
       title = "Degree Distribution by Accuracy of Predicted Interactions") +
  scale_x_discrete(labels = c("FALSE" = "No", "TRUE" = "Yes"))

# Alternatively, you might compute summary statistics:
summary_stats <- analysis_df %>%
  group_by(accurate) %>%
  summarise(mean_degree = mean(degree),
            median_degree = median(degree),
            n = n())
print(summary_stats)

## ---- new heatmap ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))

df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

# 1. Create island ID
df <- df %>%
  # filter(removed == 1) %>% mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(island_id = paste(train_layer, test_layer, sep = "_"))

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

table(df_island$observed)

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
  filter(!(avg_prop == 0 & avg_sigm_predicted < 0.5))

# 4. Order the plants (node_from) by their degree (number of interactions)
#    Here, we count how many unique pollinators each plant has in our aggregated data.
plant_degree <- df_summary %>%
  group_by(node_from) %>%
  summarise(degree = n(), .groups = "drop") %>%
  arrange(degree)

# Set the order of node_from (plants) as a factor so that the y-axis is ordered by degree.
df_summary$node_from <- factor(df_summary$node_from, levels = plant_degree$node_from)

# Optionally, you may also order pollinators (node_to) alphabetically to avoid duplicates.
# df_summary_not_all$node_to <- factor(df_summary_not_all$node_to, levels = sort(unique(df_summary_not_all$node_to)))
pollinator_degree <- df_summary %>%
  group_by(node_to) %>%
  summarise(pollinator_degree = n(), .groups = "drop") %>%
  arrange(desc(pollinator_degree))

# Set node_to factor levels based on pollinator degree (highest degree first)
df_summary$node_to <- factor(df_summary$node_to, levels = pollinator_degree$node_to)

### ---- identify always correct predictions ----
# First, check prediction correctness per iteration within each island
df_iter <- df %>%
  group_by(node_from, node_to, island_id, itr) %>%
  summarise(
    iter_correct = as.integer((original_links == 1 & sigm_predicted > 0.5) |
                                (original_links == 0 & sigm_predicted <= 0.5)),
    .groups = "drop"
  )

# For each island, the interaction is correct if every iteration was correct
df_island_t <- df_iter %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    island_correct = as.integer(all(iter_correct == 1)),
    .groups = "drop"
  )

# Finally, an interaction is always correct if it is correct in every island
df_interaction_summary <- df_island_t %>%
  group_by(node_from, node_to) %>%
  summarise(
    always_correct = as.integer(all(island_correct == 1)),
    prop_correct = mean(island_correct),  # proportion of islands correct (for information)
    n_islands = n(),
    .groups = "drop"
  )

df_summary <- df_summary %>%
  left_join(
    df_interaction_summary %>% 
      select(node_from, node_to, always_correct) %>% 
      distinct(),
    by = c("node_from", "node_to")
  )

df_summary <- df_summary %>% left_join(plant_degree,
                                       by = c("node_from")
)
df_summary$node_from <- factor(df_summary$node_from, levels = plant_degree$node_from)

# 5. Plot the heatmap and overlay circles.
#    - The background tiles (geom_tile) show the averaged proportion observed (avg_prop).
#    - We overlay circles (geom_point) for interactions where avg_sigm_predicted > 0.6.
#      Here, we differentiate by color:
#         * Use red if avg_prop == 0 (i.e. never observed empirically)
#         * Use blue if avg_prop > 0
ggplot(df_summary, aes(x = node_to, y = node_from, fill = avg_prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plant",
    fill = "Proportion\nof islands\nobserved"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 8) 
  ) + 
  # Overlay circles for high predicted probability (avg_sigm_predicted > 0.6)
  # geom_point(
  #   data = df_summary %>% filter(avg_sigm_predicted > 0.6, avg_prop == 0),
  #   aes(x = node_to, y = node_from),
  #   color = "thistle", shape = 19, size = 3, alpha = 0.6
  # ) +
  geom_point(
    data = df_summary %>% filter(always_correct == 1 & avg_prop != 0),
    aes(x = node_to, y = node_from),
    shape = 21,     # circle shape that supports fill and border
    fill = NA,      # no fill (empty circle)
    color = "salmon",  # outline color
    size = 4.5,       # adjust size as needed
    stroke = 1,    # adjust border thickness
    alpha = 0.7
  )  +
  tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

# First, add a new column "circle_type" to df_summary:
df_summary <- df_summary %>%
  mutate(circle_type = case_when(
    avg_sigm_predicted > 0.6 & avg_prop == 0 ~ "predicted, never observed",
    always_correct == 1 & avg_prop != 0 ~ "always predicted correctly",
    TRUE ~ NA_character_
  ))

# Then, modify your ggplot code:
ggplot(df_summary, aes(x = node_to, y = node_from, fill = avg_prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  # Add a combined geom_point layer for circle_type:
  geom_point(
    data = df_summary %>% filter(!is.na(circle_type)),
    aes(shape = circle_type, color = circle_type),
    size = 4.5, stroke = 1, fill = NA, alpha = 0.7
  ) +
  # Manually assign shapes and colors to circle_type levels:
  scale_shape_manual(
    name = "Prediction type",
    values = c("predicted, never observed" = 20, "always predicted correctly" = 21)
  ) +
  scale_color_manual(
    name = "Prediction type",
    values = c("predicted, never observed" = "thistle", "always predicted correctly" = "salmon")
  ) +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plant",
    fill = "Proportion\nof islands\nobserved"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4)
    #axis.text.y = element_text(size = 8)
  ) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) + tme

### ---- two scales ----
library(ggnewscale)

ggplot(df_summary, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop)) +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Proportion\nof islands\nobserved") +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = df_summary %>% filter(avg_prop == 0, avg_sigm_predicted > 0.6),
    aes(fill = avg_sigm_predicted),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2", 
                      name = "Average \npredicted \nprobability") +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

### ---- always predicted correctly separately ----
ggplot(df_summary, aes(x = node_to, y = node_from, fill = avg_prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plant",
    fill = "Proportion\nof islands\nobserved"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 8) 
  ) + 

  geom_point(
    data = df_summary %>% filter(always_correct == 1 & avg_prop != 0),
    aes(x = node_to, y = node_from),
    shape = 21,     # circle shape that supports fill and border
    fill = NA,      # no fill (empty circle)
    color = "palegreen3",  # outline color
    size = 4.5,       # adjust size as needed
    stroke = 1,    # adjust border thickness
    alpha = 0.7
  )  +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal"
  ) +
  tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

df_summary <- df_summary %>%
  mutate(circle = case_when(
    always_correct == 1 & avg_prop != 0 ~ "always predicted correctly",
    TRUE ~ NA_character_
  ))

# Then, modify your ggplot code:
ggplot(df_summary, aes(x = node_to, y = node_from, fill = avg_prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  # Add a combined geom_point layer for circle_type:
  geom_point(
    data = df_summary %>% filter(!is.na(circle)),
    aes(shape = circle, color = circle),
    size = 4.5, stroke = 1, fill = NA, alpha = 0.7
  ) +
  # Manually assign shapes and colors to circle_type levels:
  scale_shape_manual(
    name = "Prediction type",
    values = c("always predicted correctly" = 21)
  ) +
  scale_color_manual(
    name = "Prediction type",
    values = c("always predicted correctly" = "palegreen3")
  ) +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plant",
    fill = "Proportion\nof islands\nobserved"
  ) +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal"
  ) +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) + tme

## --- correlation with degree ----
# this should be corrected because the degree should be calculated based on all interactions (not just removed ones)
df_summary <- df_summary %>% left_join(df_interaction_summary,
                                       by = c("node_from", "node_to"))
view(df_summary)
ggplot(df_summary, aes(x = degree, y = avg_sigm_predicted)) +
  geom_point()

# Group by pollinator (node_to) and calculate the average proportion correct and average degree
df_pollinator_summary <- df_summary %>%
  group_by(node_to) %>%
  summarise(
    avg_prop_correct = mean(prop_correct, na.rm = TRUE),
    avg_degree = mean(degree, na.rm = TRUE),
    avg_sigm_predicted = mean(avg_sigm_predicted, na.rm = TRUE),
    n = n()  # optional, for diagnostic purposes
  ) %>%
  ungroup()

# Plot the average proportion correct vs. the average degree for each pollinator species
pollinator_degree1 <- ggplot(df_pollinator_summary, aes(x = avg_degree, y = avg_prop_correct)) +
  geom_point(alpha = 0.7, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Proportion correct",
    title = "Pollinators"
  ) +
  theme_minimal() + tme

# correlation
correlation <- cor.test(df_pollinator_summary$avg_sigm_predicted, df_pollinator_summary$avg_degree, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

pollinator_degree2 <- ggplot(df_pollinator_summary, aes(x = avg_degree, y = avg_sigm_predicted)) +
  geom_point(alpha = 0.7, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Average predicted probability (non-observed)",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 150, y = 0.59,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black")
# plants

# Group by pollinator (node_to) and calculate the average proportion correct and average degree
df_plant_summary <- df_summary %>%
  group_by(node_from) %>%
  summarise(
    avg_prop_correct = mean(prop_correct, na.rm = TRUE),
    avg_degree = mean(degree, na.rm = TRUE),
    n = n()  # optional, for diagnostic purposes
  ) %>%
  ungroup()

df_never_observed <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_plant_summary <- df_never_observed %>%
  left_join(df_plant_summary, by = "node_from")

# Plot the average proportion correct vs. the average degree for each plant species

# correlation
correlation_plants <- cor.test(df_plant_summary$count_never_observed, df_plant_summary$avg_degree, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 3)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree1 <- ggplot(df_plant_summary, aes(x = avg_degree, y = count_never_observed)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Number of predicted, non-observed interactions",
    title = "Plants"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 90, y = 165,   # Adjust depending on your data range
           label = label_text_plants,
           size = 5,
           color = "black")

# for pollinators
df_never_observed_poll <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_pollinator_summary <- df_never_observed_poll %>%
  left_join(df_pollinator_summary, by = "node_to")

# Plot the average proportion correct vs. the average degree for each plant species

# correlation
correlation_poll <- cor.test(df_pollinator_summary$count_never_observed, df_pollinator_summary$avg_degree, use = "complete.obs", method = "pearson")
correlation_poll
# Extract correlation coefficient and p-value
r_value_poll <- round(correlation_poll$estimate, 3)
p_value_poll <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value_poll, ", p = ", p_value_poll)

poll_degree <- ggplot(df_pollinator_summary, aes(x = avg_degree, y = count_never_observed)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Number of predicted, non-observed interactions",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 130, y = 26,   # Adjust depending on your data range
           label = label_text,
           size = 5,
           color = "black")

final_plot <- combine_plots(plant_degree, poll_degree)

# checking correlation with degree
library(dplyr)

# Step 1: Filter the data
df_filtered <- df %>%
  filter(itr == 1, original_links == 1)

# Step 2a: Calculate degree for each plant species (node_from)
plant_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(degree = n(), .groups = "drop")

# Average plant degree by train_layer and test_layer
avg_plant_degree <- plant_degree %>%
  group_by(node_from) %>%
  summarise(avg_plant_degree = mean(degree), .groups = "drop")

# Step 2b: Calculate degree for each pollinator species (node_to)
pollinator_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(degree = n(), .groups = "drop")

# Average pollinator degree by train_layer and test_layer
avg_pollinator_degree <- pollinator_degree %>%
  group_by(node_to) %>%
  summarise(avg_pollinator_degree = mean(degree), .groups = "drop")

# Optional: Merge the two average degree data frames into one
avg_degrees <- full_join(avg_plant_degree, avg_pollinator_degree, 
                         by = c("train_layer", "test_layer"))

# Print the results
print(avg_degrees)

## ---- try again ----
df_summary2 <- df_summary %>% left_join(pollinator_degree, by = "node_to")
df_summary2 <- df_summary %>% left_join(plant_degree, by = "node_from")

df_never_observed_poll <- df_summary2 %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_summary2 <- df_summary2 %>%
  left_join(df_never_observed_poll, by = "node_to")

# Plot the average proportion correct vs. the average degree for each plant species

# correlation
correlation <- cor.test(df_summary2$count_never_observed, df_summary2$pollinator_degree, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_summary2, aes(x = pollinator_degree, y = count_never_observed)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Number of predicted, non-observed interactions",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 130, y = 26,   # Adjust depending on your data range
           label = label_text,
           size = 5,
           color = "black")

## plants

df_never_observed_plants <- df_summary2 %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_from) %>%
  summarise(count_never_observed_plants = n(), .groups = "drop")

df_summary2 <- df_summary2 %>%
  left_join(df_never_observed_plants, by = "node_from")

# Plot the average proportion correct vs. the average degree for each plant species

# correlation
correlation <- cor.test(df_summary2$count_never_observed_plants, df_summary2$degree, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_summary2, aes(x = degree, y = count_never_observed_plants)) +
  geom_point(alpha = 0.6, size = 2, color = "palegreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Degree",
    y = "Number of predicted, non-observed interactions",
    title = "Plants"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 130, y = 26,   # Adjust depending on your data range
           label = label_text,
           size = 5,
           color = "black")



