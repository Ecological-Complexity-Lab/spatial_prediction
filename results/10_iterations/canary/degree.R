library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(gridExtra)

source("~/Documents/github/softimpute/results/useful_for_plotting.R")


### ---- correlation between never observed interactions and degree ----
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

df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

## ---- degree ----
# Step 1: Filter the data
df_filtered <- df %>%
  filter(itr == 1, original_links == 1)

# Step 2a: Calculate degree for each plant species (node_from)
plant_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(plant_degree = n(), .groups = "drop")

# Average plant degree by train_layer and test_layer
avg_plant_degree <- plant_degree %>%
  group_by(node_from) %>%
  summarise(avg_plant_degree = mean(plant_degree), .groups = "drop")

# try with overall
overall_plant_degree <- df_filtered %>% 
  group_by(node_from) %>% 
  summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")

overall_poll_degree <- df_filtered %>% 
  group_by(node_to) %>% 
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

# Step 2b: Calculate degree for each pollinator species (node_to)
pollinator_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(poll_degree = n(), .groups = "drop")

# Average pollinator degree by train_layer and test_layer
avg_pollinator_degree <- pollinator_degree %>%
  group_by(node_to) %>%
  summarise(avg_pollinator_degree = mean(poll_degree), .groups = "drop")

# # # Optional: Merge the two average degree data frames into one
# avg_degrees <- full_join(avg_plant_degree, avg_pollinator_degree,
#                          by = c("train_layer", "test_layer"))
# 
# # Print the results
# print(avg_degrees)

# 1. Create island ID
df <- df %>%
  # filter(removed == 1) %>% mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(sigm_predicted = sigmoid(predicted_values)) %>% 
  mutate(island_id = paste(train_layer, test_layer, sep = "_"))

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

df_never_observed <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_never_observed_poll <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

# # degree calculation 
# plant_degree <- df_summary %>%
#   group_by(node_from) %>%
#   summarise(degree = n(), .groups = "drop") %>%
#   arrange(degree)
# 
# pollinator_degree <- df_summary %>%
#   group_by(node_to) %>%
#   summarise(pollinator_degree = n(), .groups = "drop") %>%
#   arrange(desc(pollinator_degree))

# df_plant_summary <- df_summary %>%
#   group_by(node_from) %>%
#   summarise(
#     avg_degree = mean(degree, na.rm = TRUE),
#     n = n()  # optional, for diagnostic purposes
#   ) %>%
#   ungroup()
# 
# df_summary <- df_never_observed %>%
#   left_join(df_summary, by = "node_from")

never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")

# for overall degree

never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")

never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")

df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

df_to_correlate <- never_poll_degree_overall
df_to_correlate$x <- df_to_correlate$overall_poll_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

# correlation
correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 3)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree1 <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, non-observed interactions",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = 19, y = 27,   # Adjust depending on your data range
           label = label_text_plants,
           size = 5,
           color = "black")
plant_degree1

# for pollinators
df_never_observed_poll <- df_summary %>%
  filter(avg_prop == 0, avg_sigm_predicted > 0.5) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_summary <- df_never_observed %>%
  left_join(df_summary, by = "node_to")

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
