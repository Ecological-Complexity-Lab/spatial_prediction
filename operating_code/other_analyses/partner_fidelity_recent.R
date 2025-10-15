### ---- partner fidelity correlation with evaluators ----
# this part of the code was removed from Abramov_et_al_spatial_prediction_analysis.R

# filter existing interactions in cases where train = test layer (we want to calculate the average similarity in partner composition of each species in different locations)
df_fidelity <- df %>% filter(train_layer == test_layer) %>% 
  filter(original_links != 0) %>% filter(itr == 1) # all iterations are the same for existing links that were not removed

#### ---- plants fidelity ----
# for each plant (node_from) and layer, gather the pollinators (node_to).
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
      # If the plant is only in one layer, there are no pairs, so return NA. select only species the occur in at least 3 islands
      if (n_layers < 3) {
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

# we filter out single-island predictions since what is relevant here is prediction across space
df_off <- df %>% 
  filter(train_layer != test_layer) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0)) %>% 
  mutate(original_binary = if_else(original_links > 0, 1, 0))

df_merged <- df_off %>%
  left_join(df_sorensen, by = "node_from")

#### ---- pollinators fidelity ----
df_pollinators <- df_fidelity %>%
  distinct(node_to, train_layer, node_from) %>%
  group_by(node_to, train_layer) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# 2. For each pollinator, compute mean Sorensen similarity across all pairs of layers.
df_sorensen_pollinators <- df_pollinators %>%
  group_by(node_to) %>%
  summarise(
    mean_sorensen_pollinators = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 3) {
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
  left_join(df_sorensen_pollinators, by = "node_to") # joint plant and pollinator fidelities

#### ---- per-species f1 ----
# how many removed 1s and 0s do we have for each species?
# For plants:
removed_count_plants <- df_merged %>%
  filter(removed == 1) %>% 
  group_by(plant = node_from) %>%
  summarise(
    n_removed_0    = sum(original_binary == 0, na.rm = TRUE),
    n_removed_1    = sum(original_binary == 1, na.rm = TRUE),
  ) %>%
  arrange(plant)

# For pollinators:
removed_count_polls <- df_merged %>%
  filter(removed == 1) %>% 
  group_by(pollinator = node_to) %>%
  summarise(
    n_removed_0    = sum(original_binary == 0, na.rm = TRUE),
    n_removed_1    = sum(original_binary == 1, na.rm = TRUE),
  ) %>%
  arrange(pollinator)

# calculate balanced f1 for each species by sampling its removed links 50 times in a balanced subset and calculating f1 for each
df_clean <- df_merged %>%
  filter(removed == 1) %>% 
  filter(!is.na(original_binary), !is.na(predicted_bin_sigm))

set.seed(42)

# 3) For pollinators (node_to), bootstrap 50× per species
df_balanced_f1_polls <- df_clean %>%
  group_by(species = node_to) %>%
  group_modify(~ {
    f1s <- replicate(50, compute_balanced_f1(.x))
    tibble(
      mean_f1 = mean(f1s, na.rm = TRUE),
      sd_f1   = sd(  f1s, na.rm = TRUE)
    )
  }) %>%
  ungroup()

# 4) Do the same for plants (node_from)
df_balanced_f1_plants <- df_clean %>%
  group_by(species = node_from) %>%
  group_modify(~ {
    f1s <- replicate(50, compute_balanced_f1(.x))
    tibble(
      mean_f1 = mean(f1s, na.rm = TRUE),
      sd_f1   = sd(  f1s, na.rm = TRUE)
    )
  }) %>%
  ungroup()

# join them together
df_for_plot_balanced <- bind_rows(
  df_balanced_f1_polls %>%
    left_join(df_sorensen_pollinators, by = c("species" = "node_to")) %>%
    mutate(avg_sorensen_plants = NA_real_),
  
  df_balanced_f1_plants %>%
    left_join(df_sorensen, by = c("species" = "node_from")) %>%
    mutate(avg_sorensen_pollinators = NA_real_)
)

#### ---- Fig. S13: plot fidelity ----
balanced_f1_fidelity <- make_full_correlation_plot(
  data            = df_for_plot_balanced,
  evaluator       = "mean_f1",
  pollinator_x    = "mean_sorensen_pollinators",
  plant_x         = "mean_sorensen_plants",
  shared_x_lab    = "Mean partner fidelity (Sørensen)",
  shared_y_lab    = "F1 score per species"
)

# pdf("balanced_f1_fidelity.pdf", width = 6, height = 4)
# grid::grid.draw(balanced_f1_fidelity)
# dev.off()

