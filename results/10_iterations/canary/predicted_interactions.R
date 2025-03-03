library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

df <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')

# 1. Keep only rows where train_layer == test_layer (i.e., “islands”)
df_filtered <- df %>%
  filter(train_layer == test_layer)

# if you want to work with all combinations of train and test, switch df_filtered with df

# 2. Create an island identifier
df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_"))

df <- df %>% filter(removed == 1)

# 3. For each interaction, compute the proportion of islands in which it was observed
df_prop <- df %>%
  group_by(node_from, node_to) %>%
  summarize(
    prop_islands_observed = mean(original_links == 1, na.rm = TRUE),
    # This can be handy to see how many islands contributed:
    n_islands = n_distinct(island_id)
  ) %>%
  ungroup()

# 4. (Optional) If you only want interactions that are NOT observed in all islands:
df_prop_not_all <- df_prop %>%
  filter(prop_islands_observed < 1)

# 5. Plot as a heatmap
#    - If you want a complete matrix, you can pivot_wider() to fill missing combos
#      or you can directly plot the long format:
ggplot(df_prop, aes(x = node_to, y = node_from, fill = prop_islands_observed)) +
  geom_tile() +
  # You can tweak the color scale to suit your taste:
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    x = "Plant",
    y = "Pollinator",
    fill = "Prop. observed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x labels if they're long
  )

# add predicted probabilities
df <- df %>% mutate(sigm_predicted = sigmoid(predicted_values))

ggplot(df, aes(x = node_to, y = node_from, fill = sigm_predicted)) +
  geom_tile() +
  # You can tweak the color scale to suit your taste:
  scale_fill_gradient(low = "white", high = "thistle3") +
  theme_minimal() +
  labs(
    x = "Plant",
    y = "Pollinator",
    fill = "Predicted probability"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # rotate x labels if they're long
  ) + tme

# try to plot both

# 1. Merge the two data frames by (node_from, node_to)
#    Suppose df_prop has columns: node_from, node_to, prop_islands_observed
#    and df_pred has columns: node_from, node_to, predicted_probability

df_combined <- df_prop %>%
  left_join(df, by = c("node_from", "node_to"))

# 2. Plot
ggplot(df_combined, aes(x = node_to, y = node_from)) +
  # a) Heatmap background for proportion observed
  geom_tile(aes(fill = prop_islands_observed)) +
  
  # b) Points for predicted prob > 0.5, color-coded by observed proportion
  #    RED if observed proportion == 0, BLUE if observed proportion == 1
  geom_point(
    data = df_combined %>% filter(sigm_predicted > 0.6, prop_islands_observed == 0),
    color = "salmon", shape = 19, size = 3, alpha = 0.2
  ) +
  geom_point(
    data = df_combined %>% filter(sigm_predicted > 0.6, prop_islands_observed > 0),
    color = "thistle", shape = 19, size = 3, alpha = 0.1
  ) +

  # Adjust fill scale, theme, etc.
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    x = "Pollinator",
    y = "Plants",
    fill = "Prop. Observed"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(igraph)
library(dplyr)

# Assume df_combined has: node_from, node_to, prop_islands_observed, predicted_probability
# Create an edge list:
edges <- df_combined %>% 
  select(node_from, node_to, prop_islands_observed, sigm_predicted)

# Create an igraph object (undirected or directed based on your context)
g <- graph_from_data_frame(edges, directed = FALSE)

# Assign colors to edges based on your conditions:
edge_colors <- ifelse(
  edges$prop_islands_observed == 0 & edges$sigm_predicted > 0.5, "salmon",
  ifelse(edges$prop_islands_observed > 0 & edges$sigm_predicted > 0.5, "steelblue", "grey")
)
E(g)$color <- edge_colors

# Optionally, you can add the predicted_probability as an edge attribute:
E(g)$sigm_predicted <- edges$sigm_predicted

# Choose a layout (e.g., layout_with_fr for a force-directed layout)
layout <- layout_with_fr(g)

# Plot the network
plot(g, layout = layout,
     vertex.label = V(g)$name,
     vertex.size = 1,  # adjust size as needed
     main = "Network of Interactions with Highlighted Edges")

# Load required packages
library(bipartite)
library(dplyr)

# Suppose your combined data frame is called df_combined and looks like:
#   node_from                node_to                prop_islands_observed   predicted_probability
# 1 "Thomisus_cf_onustus"    "Euphorbia_balsamifera_m"   1                    0.000459518
# 2 "Thysanoptera_sp_2"      "Phagnalon_saxatile"        0                    0.75
# ... (your data)

# Filter to keep interactions with predicted_probability > 0.5 in ALL iterations and layer combos,
# then summarize to get average predicted probability and average proportion observed.
df_summary <- df_combined %>%
  group_by(node_from, node_to) %>%
  # Keep groups where every observation meets the condition
  filter(all(sigm_predicted > 0.5)) %>%
  summarise(
    avg_predicted_probability = mean(sigm_predicted),
    avg_prop_islands_observed = mean(prop_islands_observed),
    n = n(),  # number of rows contributing, useful for debugging or further filtering
    .groups = "drop"
  )

# Now df_summary can be used in your bipartite graph.
head(df_summary)
# 1. Create the incidence matrix.
#    Rows: node_from, Columns: node_to.
inc_mat <- xtabs(~ node_from + node_to, data = df_summary)

# 2. Create a matrix to hold the colors for each interaction.
#    Initialize all cells to a default color (e.g., "grey").
color_mat <- matrix("grey", 
                    nrow = nrow(inc_mat), 
                    ncol = ncol(inc_mat), 
                    dimnames = dimnames(inc_mat))

# 3. Loop through each row in df_combined to update the color matrix based on conditions.
for(i in seq_len(nrow(df_summary))) {
  current <- df_combined[i, ]
  # Convert node names to characters (if factors)
  from <- as.character(current$node_from)
  to <- as.character(current$node_to)
  
  if(current$sigm_predicted > 0.5) {
    if(current$prop_islands_observed == 0) {
      color_mat[from, to] <- "salmon"
    } else if(current$prop_islands_observed > 0) {
      color_mat[from, to] <- "steelblue"
    }
  }
}

# 4. Plot the bipartite network using plotweb from the bipartite package.
plotweb(inc_mat,
        col.interaction = color_mat, # Use our custom color matrix for links
        labsize = 1,                 # Adjust label size as needed
        text.rot = 45) #,               # Rotate text for better readability


visweb(inc_mat)
