df_merged_rem <- df_merged %>% filter(removed == 1) %>% 
  mutate(predicted_sigm = sigmoid(predicted_values))

library(dplyr)

library(dplyr)

# Filter only rows where there is an interaction (original_links == 1)
df_itr <- df %>% filter(itr == 1)
edges <- df_itr %>% 
  filter(original_links == 1) 

# --------------------------
# 1. Calculate Total Degree
# --------------------------
# Create an edge list for both directions: 
# from node_from to node_to and from node_to to node_from.
edges_from <- edges %>% 
  select(species = node_from, partner = node_to)
edges_to <- edges %>% 
  select(species = node_to, partner = node_from)

# Combine the two to capture all interactions
all_edges <- bind_rows(edges_from, edges_to) %>% 
  distinct()  # ensuring uniqueness

# For each species, count the unique partners
total_degree <- all_edges %>% 
  group_by(species) %>% 
  summarize(total_degree = n_distinct(partner)) %>% 
  ungroup()

# ---------------------------
# 2. Calculate Average Degree
# ---------------------------
# For each layer combination, first compute the degree per species
# (again, considering both directions).
# We'll get an edge list for each layer, then count per species.
edges_layer <- edges %>%
  group_by(train_layer, test_layer) %>%
  do({
    layer_edges <- .
    edges_from <- layer_edges %>% select(species = node_from, partner = node_to)
    edges_to <- layer_edges %>% select(species = node_to, partner = node_from)
    bind_rows(edges_from, edges_to) %>% distinct()
  }) %>% ungroup()

# Now, for each layer combination and species, calculate degree
degree_per_layer <- edges_layer %>%
  group_by(train_layer, test_layer, species) %>%
  summarize(degree = n_distinct(partner)) %>%
  ungroup()

# Finally, average the degree for each species across layers
avg_degree <- degree_per_layer %>%
  group_by(species) %>%
  summarize(avg_degree = mean(degree)) %>%
  ungroup()

# Optionally, merge the two results to have a summary table per species:
degree_summary <- total_degree %>%
  left_join(avg_degree, by = "species")

# Inspect the result
print(degree_summary)

library(dplyr)

df_joined <- df %>%
  left_join(degree_summary, by = c("node_to" = "species")) %>%
  left_join(degree_summary, by = c("node_from" = "species"), suffix = c("_to", "_from"))

df_joined <- df_joined %>% mutate(predicted_sigm = sigmoid(predicted_values))

# count interactions that were never observed yet consistantly predicted as existing

library(dplyr)

# Identify interactions that always meet the criteria
valid_interactions <- df_joined %>%
  group_by(node_to, node_from) %>%
  # Check if every row in the group satisfies the condition:
  # note that using all() inside summarize ensures that the condition holds in every iteration.
  summarize(
    meets_condition = all(original_links == 0 & predicted_sigm > 0.5),
    n_layers = n(),  # you can check how many rows contributed per interaction
    .groups = "drop"
  ) %>%
  # Filter to keep only interactions where the condition was met in every instance
  filter(meets_condition)

# Count how many unique interactions meet the criteria
interaction_count <- nrow(valid_interactions)
print(interaction_count)


valid_links <- df_joined %>% filter(original_links == 0 & predicted_sigm > 0.5)
valid_links_all <- df_joined %>% filter(all(original_links == 0 & predicted_sigm > 0.5))
