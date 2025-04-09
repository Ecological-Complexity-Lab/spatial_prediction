# ---- partner fidelity ----
# here we calculate the "partner fidelity" (measured by Sorensen similarity between the partner composition of each plant and pollinator) and compare it to the accuracy of predictions

# themes
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load libraries ----
library(dplyr)
library(tidyverse)

## ---- load data ----
df <- read_csv('nonbinary_equal_0_1_removal_60_1.csv') # canary islands. replace with the nonbinary version or just add the weights from it
df <- df %>% filter(k == 2 & lambda == 0.1)
#df <- read_csv('test_df_itr_1_60_binary.csv') # test on 1 iteration

# filter out cases in which train = test layer
df_fidelity <- df %>% filter(train_layer == test_layer) %>% 
  filter(original_links == 1) %>% filter(itr == 1)

## ---- plants ----
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

# for pollinators

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

write.csv(df_sorensen_pollinators, "sorensen_pollinators_site.csv")
write.csv(df_sorensen, "sorensen_plants_site.csv")

## ---- add distance between predicted and observed values ----
df_merged <- df_merged %>% mutate(predicted_value_sigm = sigmoid(predicted_values))
df_merged <- df_merged %>% mutate(predicted_sigm_obs = predicted_value_sigm - original_links)
df_merged <- df_merged %>% mutate(pred_obs_abs = abs(predicted_sigm_obs))

ggplot(df_merged, aes(x = mean_sorensen_plants, y = pred_obs_abs)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "Predicted - observed",
       title = "Sigmoid predicted value vs. mean partner fidelity of pollinators (Sorensen)") +
  tme
