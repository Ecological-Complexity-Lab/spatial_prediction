### ---- degree impact on link assignment ----
# here we examine if the algorithm assigns more links to species with higher degree.

#### ---- calculate overall degree per species ----
# filter only existing links
df_filtered <- df %>%
  filter(itr == 1, original_links != 0) # filter existing interactions

# calculate degree for each plant species (node_from)
plant_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(plant_degree = n(), .groups = "drop")

# # average plant degree by train_layer and test_layer
# avg_plant_degree <- plant_degree %>%
#   group_by(node_from) %>%
#   summarise(avg_plant_degree = mean(plant_degree), .groups = "drop")

overall_plant_degree <- df_filtered %>% 
  group_by(node_from) %>% 
  summarise(overall_plant_degree = length(unique(node_to)), .groups = "drop")

# calculate degree for each pollinator species (node_to)
pollinator_degree <- df_filtered %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(poll_degree = n(), .groups = "drop")

# # average pollinator degree by train_layer and test_layer
# avg_pollinator_degree <- pollinator_degree %>%
#   group_by(node_to) %>%
#   summarise(avg_pollinator_degree = mean(poll_degree), .groups = "drop")

overall_poll_degree <- df_filtered %>% 
  group_by(node_to) %>% 
  summarise(overall_poll_degree = length(unique(node_from)), .groups = "drop")

#### ---- Fig. 3c: plot degree vs. number of never observed interactions ----
df_test <- df %>%
  mutate(island_id = paste(test_layer)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

all(df_test$test_layer == as.numeric(df_test$island_id)) # check island ids

# for each island and interaction, determine if the interaction was observed.
# we use `any(original_links == 1)` so that if the interaction is observed in at least one iteration, we count it.
df_island <- df_test %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )

table(df_island$observed)



# now, for each unique interaction, compute:
# - The proportion of islands where it was observed.
# - The average predicted probability (averaged over islands).
# df_summary <- df_island %>%
#   group_by(node_from, node_to) %>%
#   summarise(
#     avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
#     avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
#     n_islands = n(),  # number of islands contributing
#     .groups = "drop"
#   )

df_summary <- df_island %>%
  group_by(node_from, node_to) %>%
  summarise(
    n_islands_observed = sum(observed == 1L, na.rm = TRUE),
    avg_prop_isl       = n_islands_observed / 7,   # fixed denominator
    
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    .groups = "drop"
  ) # add a column to count islands for which the interaction was predicted

df_summary <- df_summary %>%
  # join pollinator degree by node_to
  left_join(overall_poll_degree, by = "node_to") %>%
  # join plant degree by node_from
  left_join(overall_plant_degree, by = "node_from")


view(df_summary)

# stopped here. fix the way island_sigm_predicted is calculated to include only observed == 0? but when we predicted we didn't care about all of the other islands. so it does not make sense.

# # filter the interactions that were never observed throughout the data set
# df_never_observed <- df_summary %>%
#   filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
#   group_by(node_from) %>%
#   summarise(count_never_observed = n(), .groups = "drop") #plants
# 
# df_never_observed_poll <- df_summary %>%
#   filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
#   group_by(node_to) %>%
#   summarise(count_never_observed = n(), .groups = "drop")
# 
# never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")
# 
# # for relating to overall degree
# 
# never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")
# 
# never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")
# 
# # calculate correlation
# df_to_correlate <- never_plants_degree_overall
# df_to_correlate$x <- df_to_correlate$overall_plant_degree
# df_to_correlate$y <- df_to_correlate$count_never_observed
# 
# correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
# correlation_plants
# 
# # Extract correlation coefficient and p-value
# r_value <- round(correlation_plants$estimate, 2)
# p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
# label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)
# 
# # specify species you want to label, if any
# from_label <- c("Euphorbia_balsamifera_m",
#                 "Euphorbia_balsamifera_f",
#                 "Launaea_arborescens")
# 
# # build a little data‐frame just for the labels
# df_labels <- df_to_correlate %>%
#   filter(node_from %in% from_label) %>%
#   distinct(node_from, .keep_all = TRUE) %>%
#   mutate(
#     label = paste0(
#       'italic("',
#       gsub("_", " ", node_from),
#       '")'
#     )
#   )
# 
# # plot
# 
# plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
#   geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
#   geom_smooth(method = "lm", se = FALSE, color = "navy") +
#   labs(
#     x = "Overall degree",
#     y = "Number of predicted, \nnon-observed interactions",
#     title = paste("Plants:", label_text_plants)   # <--- add label in title
#   ) +
#   theme_minimal() + tme
# 
# # add repel‐text layer
# plant_degree <- plant_degree +
#   geom_text_repel(
#     data      = df_labels,
#     aes(label = label),
#     parse     = TRUE,       # interpret the label as an expression
#     size      = 4,          # tweak text size as needed
#     box.padding   = 0.35,   # how much to push labels away from each other
#     point.padding = 0.5,    # how much to push labels away from the points
#     nudge_y       = -0.2     # optional small shift upward
#   )
# 
# # same for pollinators
# df_to_correlate_poll <- never_poll_degree_overall
# df_to_correlate_poll$x <- df_to_correlate_poll$overall_poll_degree
# df_to_correlate_poll$y <- df_to_correlate_poll$count_never_observed
# 
# correlation_poll <- cor.test(df_to_correlate_poll$x, df_to_correlate_poll$y, use = "complete.obs", method = "pearson")
# correlation_poll
# 
# # Extract correlation coefficient and p-value
# r_value <- round(correlation_poll$estimate, 2)
# p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
# label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)
# 
# poll_degree <- ggplot(df_to_correlate_poll, aes(x = x, y = y)) +
#   geom_point(alpha = 0.6, size = 2, color = "rosybrown2") +
#   geom_smooth(method = "lm", se = FALSE, color = "navy") +
#   labs(
#     x = "Overall degree",
#     y = "Number of predicted, \nnon-observed interactions",
#     title = paste("Pollinators:", label_text_polls)   # <--- add label in title
#   ) +
#   theme_minimal() + tme
# 
# # # Create the figure
# # final_plot <- combine_plots(plant_degree, poll_degree) # Fig. 3c
# 
# # Your two plots (remove individual axis labels)
# plant_degree_clean <- plant_degree +
#   labs(x = NULL, y = NULL)
# 
# poll_degree_clean <- poll_degree +
#   labs(x = NULL, y = NULL)
# 
# # Add bottom label with padding
# final_plot <- plot_grid(
#   # main plots
#   plot_grid(plant_degree_clean, poll_degree_clean, ncol = 2, align = "hv"),
#   # x label
#   ggdraw() + draw_label("Overall degree", fontface = "bold", size = 16),
#   ncol = 1,
#   rel_heights = c(1, 0.08)  # second element is space for x-axis label
# )
# 
# # Add y label with padding
# final_plot <- plot_grid(
#   ggdraw() + draw_label("Number of predicted,\nnon-observed links",
#                         angle = 90, fontface = "bold", size = 16),
#   final_plot,
#   ncol = 2,
#   rel_widths = c(0.08, 1)   # first element is space for y-axis label
# )
# 
# final_plot # fig. 3c

# # Save to PDF
# pdf("degree_unobserved_links.pdf", width = 10, height = 7)  # adjust size as needed
# grid::grid.draw(final_plot)
# dev.off()

# png(
#   filename = "degree_unobserved_links.png",
#   width    = 7,
#   height   = 4,
#   units    = "in",
#   res      = 300,
#   family   = "Helvetica"
# )
# 
# # --- your plotting code here ---
# print(final_plot)
# 
# dev.off()


### ---- Fig. 3a: mapping never-observed links ----
# here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of cases in which these interactions were observed.

# # order species by their degree
# df_summary <- df_summary %>% left_join(overall_poll_degree, by="node_to")
# df_summary <- df_summary %>% left_join(overall_plant_degree, by="node_from")
# 
# # determine the order of species in the plot based on their degree
plant_order <- df_summary %>%
  distinct(node_from, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(node_from)

poll_order <- df_summary %>%
  distinct(node_to, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(node_to)
# 
# # calculate proportion of islands in which each interaction occurs, rather than island pairs
# df_filtered_self <- df_filtered %>% filter(train_layer == test_layer)
# 
# result_count <- df_filtered_self %>%
#   # group by the interaction
#   group_by(node_to, node_from) %>%
#   summarise(
#     # count unique test layers for this interaction
#     n_test_layers = n_distinct(test_layer),
#     .groups = "drop"
#   ) %>%
#   # calculate the proportion
#   mutate(
#     prop_test_layers = n_test_layers / n_distinct(df$test_layer)
#   )
# 
# print(result_count)
# 
# # join to df_summary
# 
# df_summary <- df_summary %>%
#   left_join(result_count %>% select(node_to, node_from, prop_test_layers),
#             by = c("node_to", "node_from")) %>%
#   # replace avg_prop with the calculated proportion
#   mutate(avg_prop_isl = if_else(is.na(prop_test_layers), 0, prop_test_layers)) %>%
#   select(-prop_test_layers)  # remove helper column if not needed
# 
# all((df_summary$avg_prop == 0) == (df_summary$avg_prop_isl == 0)) # check

# reset the levels for the species factors
df_summary$node_from <- factor(df_summary$node_from, levels = plant_order)
df_summary$node_to   <- factor(df_summary$node_to, levels = poll_order)

map_missing_links <- ggplot(df_summary, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop_isl)) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = df_summary %>% filter(avg_prop_isl == 0, avg_sigm_predicted > best_discrete_threshold),
    aes(fill = avg_sigm_predicted),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_missing_links)
# 
# pdf(
#   file   = "map_missing_links.pdf",
#   width  = 11,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(map_missing_links)
# dev.off()     # close the file

# png(
#   filename = "map_missing_links.png",
#   width    = 12,
#   height   = 7,
#   units    = "in",
#   res      = 300,
#   family   = "Helvetica"
# )
# 
# # --- your plotting code here ---
# print(map_missing_links)
# 
# dev.off()

map_missing_links_all_noleg <- map_missing_links +
  theme(legend.position = "none")


map_missing_links_blue_noleg <- map_missing_links +
  theme(legend.position = "none")

# png(
#   filename = "map_missing_links_blue.png",
#   width    = 12,
#   height   = 6,
#   units    = "in",
#   res      = 300,
#   family   = "Helvetica"
# )
# 
# # --- your plotting code here ---
# print(map_missing_links_blue_noleg)
# 
# dev.off()

### ---- Fig. S11: difference in links predicted with/without external data ----
# this analysis shows us which links (and how many) were predicted only using external data, single-island data or combination of both.
# df_island_sep <- df_island %>%

# only within-island data

df_test_diag <- df_test %>%
  filter(island_id == train_layer)

df_island_diag <- df_test_diag %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )


df_summary_diag <- df_island_diag %>%
  group_by(node_from, node_to) %>%
  summarise(
    n_islands_observed = sum(observed == 1L, na.rm = TRUE),
    avg_prop_isl       = n_islands_observed / 7,   # fixed denominator
    
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_diag <- df_summary_diag %>%
  # join pollinator degree by node_to
  left_join(overall_poll_degree, by = "node_to") %>%
  # join plant degree by node_from
  left_join(overall_plant_degree, by = "node_from")


view(df_summary_diag)

# df_summary_diag <- df_island_diag %>%
#   group_by(node_from, node_to) %>%
#   summarise(
#     avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
#     avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
#     n_islands = n(),  # number of islands contributing
#     .groups = "drop"
#   )

# only with extrnal data

df_test_offs <- df_test %>%
  filter(island_id != train_layer)

df_island_offs <- df_test_offs %>%
  group_by(node_from, node_to, island_id) %>%
  summarise(
    observed = as.integer(any(original_links != 0)),
    # For predicted_prob_sigm, you might take the average across iterations per island.
    island_sigm_predicted = mean(predicted_prob_sigm, na.rm = TRUE),
    .groups = "drop"
  )


df_summary_offs <- df_island_offs %>%
  group_by(node_from, node_to) %>%
  summarise(
    n_islands_observed = sum(observed == 1L, na.rm = TRUE),
    avg_prop_isl       = n_islands_observed / 7,   # fixed denominator
    
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_offs <- df_summary_offs %>%
  # join pollinator degree by node_to
  left_join(overall_poll_degree, by = "node_to") %>%
  # join plant degree by node_from
  left_join(overall_plant_degree, by = "node_from")


view(df_summary_offs)

# df_island_offs <- df_island_sep %>%
#   filter(island1 != island2)
# 
# df_summary_offs <- df_island_offs %>%
#   group_by(node_from, node_to) %>%
#   summarise(
#     avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
#     avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
#     n_islands = n(),  # number of islands contributing
#     .groups = "drop"
#   )
# 
# identical(df_summary_offs$node_from, df_summary_diag$node_from) # check

# create a combined table
diff_df <- full_join(
  df_summary_offs,
  df_summary_diag,
  by     = c("node_from","node_to"),
  suffix = c("_offs","_diag")
) 

# these are interactions that were predicted to exist in the off-diagonals but not in the diagonal
diff_df_filtered <- diff_df %>% filter(diff_df$avg_sigm_predicted_diag < best_discrete_threshold & diff_df$avg_sigm_predicted_offs >= best_discrete_threshold)

# and these are interactions that were predicted to exist in the diagonal but not in the off-diagonals
diff_df_diags <- diff_df %>% filter(avg_sigm_predicted_offs < best_discrete_threshold & avg_sigm_predicted_diag >= best_discrete_threshold)

# order species by their degree
diff_df$node_from <- factor(diff_df$node_from, levels = plant_order)
diff_df$node_to   <- factor(diff_df$node_to, levels = poll_order)

# plot the differences
df_plot <- diff_df %>%
  filter(avg_prop_isl_diag == 0 & avg_prop_isl_offs == 0) %>% # non-observed interactions
  mutate(
    sigm_cat = case_when(
      avg_sigm_predicted_diag < best_discrete_threshold &
        avg_sigm_predicted_offs > best_discrete_threshold ~ "offs↑ only",
      
      avg_sigm_predicted_offs < best_discrete_threshold &
        avg_sigm_predicted_diag > best_discrete_threshold ~ "diag↑ only",
      
      avg_sigm_predicted_offs > best_discrete_threshold &
        avg_sigm_predicted_diag > best_discrete_threshold ~ "both↑",
      
      TRUE ~ NA_character_
    )
  )

map_missing_links_diags_offs <- ggplot(df_plot, aes(x = node_to, y = node_from)) +
  
  # allow a second fill scale
  new_scale_fill() +
  geom_tile(
    data  = filter(df_plot, !is.na(sigm_cat)),
    aes(fill = sigm_cat),
    alpha = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "offs↑ only" = "salmon",
      "diag↑ only" = "plum3",
      "both↑"       = "lightsteelblue"
    ),
    na.value = NA,
    name   = "Difference in \nprediction approach",
    labels = c(
      "offs↑ only" = "Predicted only by \nadding external location",
      "diag↑ only" = "Predicted only by \nsingle location",
      "both↑"       = "Predicted by both approaches"
    )
  ) +
  
  # tidy up
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 10),
    legend.position = "bottom",
    legend.box      = "vertical"
  ) +
  # your italic‐species labels and extra theme element
  scale_y_discrete(
    labels = function(x) lapply(strsplit(x, "_"), function(y) {
      bquote(italic(.(paste(y, collapse = " "))))
    })
  ) +
  tme

map_missing_links_diags_offs


# pdf(
#   file   = "map_missing_links_diags_offs.pdf",
#   width  = 11,    # inches
#   height = 6,
#   family = "Helvetica"   # or another installed font
# )
# print(map_missing_links_diags_offs)
# dev.off()     # close the file
# 
png(
  filename = "map_missing_links_diags_offs.png",
  width    = 12,
  height   = 7,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_missing_links_diags_offs)

dev.off()

map_missing_links_noleg <- map_missing_links +
  theme(legend.position = "none")

map_missing_links_diags_offs_noleg <- map_missing_links_diags_offs +
  theme(legend.position = "none")

png(
  filename = "map_missing_links_diags_offs_noleg.png",
  width    = 12,
  height   = 6,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_missing_links_diags_offs_noleg)

dev.off()
# 
# png(
#   filename = "map_missing_links_noleg.png",
#   width    = 12,
#   height   = 6,
#   units    = "in",
#   res      = 300,
#   family   = "Helvetica"
# )
# 
# # --- your plotting code here ---
# print(map_missing_links_noleg)
# 
# dev.off()

# how many links did each category add?
df_plot %>%
  filter(avg_prop_diag == 0) %>% 
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  )

### ---- conference missing links map ----
all((diff_df$avg_prop_isl_diag) == (diff_df$avg_prop_isl_offs)) # check

map_potential_missing_links <- ggplot(diff_df, aes(x = node_to, y = node_from)) +
  # # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "steelblue",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value with within-system data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag > best_discrete_threshold & avg_sigm_predicted_offs < best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # third layer that is "a part" of the second layer: overlay additional cells from external predictions
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold & avg_sigm_predicted_diag < best_discrete_threshold),
    aes(fill = avg_sigm_predicted_offs),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # fourth layer that overlays cells predicted by both approaches
  geom_tile(
    data = diff_df %>%
      mutate(
        diag_ok = !is.na(avg_sigm_predicted_diag) & avg_sigm_predicted_diag > thr,
        offs_ok = !is.na(avg_sigm_predicted_offs) & avg_sigm_predicted_offs > thr,
        max_pred = pmax(avg_sigm_predicted_diag, avg_sigm_predicted_offs, na.rm = TRUE)
      ) %>%
      filter(avg_prop_isl_diag == 0, diag_ok, offs_ok),
    aes(fill = max_pred),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links)

map_potential_missing_links_noleg <- map_potential_missing_links +
  theme(legend.position = "none")

png(
  filename = "map_potential_missing_links_noleg.png",
  width    = 12,
  height   = 6,
  units    = "in",
  res      = 300,
  family   = "Helvetica",
  bg = "transparent"
)

# --- your plotting code here ---
print(map_potential_missing_links_noleg)

dev.off()

# only potential missing links

map_potential_missing_links_orange <- ggplot(diff_df, aes(x = node_to, y = node_from)) +
  # # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "white",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +

  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +

  # Second layer: overlay only cells that were never observed but have high predicted value with within-system data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag > best_discrete_threshold & avg_sigm_predicted_offs < best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # third layer that is "a part" of the second layer: overlay additional cells from external predictions
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold & avg_sigm_predicted_diag < best_discrete_threshold),
    aes(fill = avg_sigm_predicted_offs),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # fourth layer that overlays cells predicted by both approaches
  geom_tile(
    data = diff_df %>%
      mutate(
        diag_ok = !is.na(avg_sigm_predicted_diag) & avg_sigm_predicted_diag > thr,
        offs_ok = !is.na(avg_sigm_predicted_offs) & avg_sigm_predicted_offs > thr,
        max_pred = pmax(avg_sigm_predicted_diag, avg_sigm_predicted_offs, na.rm = TRUE)
      ) %>%
      filter(avg_prop_isl_diag == 0, diag_ok, offs_ok),
    aes(fill = max_pred),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "tan1", high = "tomato2",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links_orange)

map_potential_missing_links_orange_noleg <- map_potential_missing_links_orange +
  theme(legend.position = "none")

png(
  filename = "map_potential_missing_links_orange_noleg.png",
  width    = 12,
  height   = 6,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_potential_missing_links_orange_noleg)

dev.off()

## ---- External vs. within system prediction for conference ----

map_potential_missing_links_diag_offs <- ggplot(diff_df, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient) - we don't show existing links
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "white",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but predicted only with external data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag < best_discrete_threshold &
                                avg_sigm_predicted_offs > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_offs),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "salmon", high = "salmon",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # third layer: overlay only cells that were never observed but predicted only with within-system data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs < best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "plum", high = "plum",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # fourth layer: overlay interactions that were predicted using both approaches
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "lightsteelblue", high = "lightsteelblue",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links_diag_offs)

map_potential_missing_links_diag_offs_noleg <- map_potential_missing_links_diag_offs +
  theme(legend.position = "none")

png(
  filename = "map_potential_missing_links_diag_offs_noleg_violet.png",
  width    = 12,
  height   = 6,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_potential_missing_links_diag_offs_noleg)

dev.off()

diff_df_salmon <- diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag < best_discrete_threshold &
                                        avg_sigm_predicted_offs > best_discrete_threshold)

unique(diff_df_salmon$node_to)

diff_df_violet <- diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs < best_discrete_threshold &
                                       avg_sigm_predicted_diag > best_discrete_threshold)


diff_df_blue <- diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold &
                                     avg_sigm_predicted_diag > best_discrete_threshold)


### ---- number of prediction instances ----
# only within-island data

df_summary_diag_count <- df_island_diag %>%
  group_by(node_from, node_to) %>%
  summarise(
    n_islands_observed = sum(observed == 1L, na.rm = TRUE),
    avg_prop_isl       = n_islands_observed / 7,   # fixed denominator
    
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands_pred = sum(island_sigm_predicted > best_discrete_threshold, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_diag_count <- df_summary_diag_count %>%
  # join pollinator degree by node_to
  left_join(overall_poll_degree, by = "node_to") %>%
  # join plant degree by node_from
  left_join(overall_plant_degree, by = "node_from")


view(df_summary_diag_count)

# only with extrnal data

df_summary_offs_counts <- df_island_offs %>%
  group_by(node_from, node_to) %>%
  summarise(
    n_islands_observed = sum(observed == 1L, na.rm = TRUE),
    avg_prop_isl       = n_islands_observed / 7,   # fixed denominator
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands_pred = sum(island_sigm_predicted > best_discrete_threshold, na.rm = TRUE),
    .groups = "drop"
  )

df_summary_offs_counts <- df_summary_offs_counts %>%
  # join pollinator degree by node_to
  left_join(overall_poll_degree, by = "node_to") %>%
  # join plant degree by node_from
  left_join(overall_plant_degree, by = "node_from")


view(df_summary_offs_counts)

# create a combined table
diff_df_count <- full_join(
  df_summary_offs_counts,
  df_summary_diag_count,
  by     = c("node_from","node_to"),
  suffix = c("_offs","_diag")
) 

# order species by their degree
diff_df_count$node_from <- factor(diff_df_count$node_from, levels = plant_order)
diff_df_count$node_to   <- factor(diff_df_count$node_to, levels = poll_order)

map_potential_missing_links_count <- ggplot(diff_df_count, aes(x = node_to, y = node_from)) +
  # # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "white",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value with within-system data
  geom_tile(
    data = diff_df_count %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag > best_discrete_threshold & avg_sigm_predicted_offs < best_discrete_threshold),
    aes(fill = n_islands_pred_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "wheat", high = "lightpink3",
                      name = "Predicted links:\number of islands\npredicted",
                      breaks = seq(0, 7, 1)) +
  
  # # Reset fill scale so the next layer can have its own gradient
  # new_scale_fill() +
  
  # third layer that is "a part" of the second layer: overlay additional cells from external predictions
  geom_tile(
    data = diff_df_count %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold & avg_sigm_predicted_diag < best_discrete_threshold),
    aes(fill = n_islands_pred_offs),
    alpha = 0.6
  ) +
  # scale_fill_gradient(low = "wheat", high = "lightpink3",
  #                     name = "Predicted links:\number of islands\npredicted",
  #                     breaks = seq(0, 7, 1)) +
  # # Reset fill scale so the next layer can have its own gradient
  # new_scale_fill() +
  
  # fourth layer that overlays cells predicted by both approaches
  geom_tile(
    data = diff_df_count %>%
      mutate(
        diag_ok = !is.na(avg_sigm_predicted_diag) & avg_sigm_predicted_diag > best_discrete_threshold,
        offs_ok = !is.na(avg_sigm_predicted_offs) & avg_sigm_predicted_offs > best_discrete_threshold,
        n_max_pred = pmax(n_islands_pred_diag, n_islands_pred_offs, na.rm = TRUE)
      ) %>%
      filter(avg_prop_isl_diag == 0, diag_ok, offs_ok),
    aes(fill = n_max_pred),
    alpha = 0.6
  ) +
  # one shared legend
  scale_fill_gradient(low = "wheat", high = "lightpink3",
                      name = "Predicted links:\nnumber of islands\npredicted",
                      breaks = seq(0, 7, 1)) +

  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links_count)

png(
  filename = "map_potential_missing_links_count.png",
  width    = 12,
  height   = 7,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_potential_missing_links_count)

dev.off()

### ---- existing predicted interactions ----
# if we want to know how each category contributed to verified existing links (that were observed in the system)
# plot the differences
df_plot_verified <- diff_df %>%
  filter(avg_prop_diag != 0) %>% # observed interactions
  mutate(
    sigm_cat = case_when(
      avg_sigm_predicted_diag < best_discrete_threshold &
        avg_sigm_predicted_offs > best_discrete_threshold ~ "offs↑ only",
      
      avg_sigm_predicted_offs < best_discrete_threshold &
        avg_sigm_predicted_diag > best_discrete_threshold ~ "diag↑ only",
      
      avg_sigm_predicted_offs > best_discrete_threshold &
        avg_sigm_predicted_diag > best_discrete_threshold ~ "both↑",
      
      TRUE ~ NA_character_
    )
  )
# # if we want to plot them
# map_existing_links_predicted <- ggplot(df_plot_verified, aes(x = node_to, y = node_from)) +
# 
#   # allow a second fill scale
#   new_scale_fill() +
#   geom_tile(
#     data  = filter(df_plot_verified, !is.na(sigm_cat)),
#     aes(fill = sigm_cat),
#     alpha = 0.6
#   ) +
#   scale_fill_manual(
#     values = c(
#       "offs↑ only" = "salmon",
#       "diag↑ only" = "plum3",
#       "both↑"       = "lightsteelblue"
#     ),
#     na.value = NA,
#     name   = "Difference in \nprediction approach",
#     labels = c(
#       "offs↑ only" = "Predicted only by \nadding external location",
#       "diag↑ only" = "Predicted only by \nsingle location",
#       "both↑"       = "Predicted by both approaches"
#     )
#   ) +
# 
#   # tidy up
#   theme_minimal() +
#   labs(x = "Pollinator", y = "Plant") +
#   theme(
#     axis.text.x  = element_blank(),
#     axis.text.y  = element_text(size = 8),
#     legend.position = "bottom",
#     legend.box      = "vertical"
#   ) +
# 
#   # your italic‐species labels and extra theme element
#   scale_y_discrete(
#     labels = function(x) lapply(strsplit(x, "_"), function(y) {
#       bquote(italic(.(paste(y, collapse = " "))))
#     })
#   ) +
#   tme
# 
# map_existing_links_predicted

# how many links did each category add?
df_plot_verified %>%
  filter(avg_prop_diag != 0) %>%
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  )

# pie chart
# 1. Count how many interactions fall into each category
df_counts <- df_plot %>%
  filter(sigm_cat != "NA") %>% 
  count(sigm_cat, name = "n") %>%
  arrange(desc(sigm_cat)) %>%          # optional: control legend/order
  mutate(
    frac = n / sum(n),                 # fraction of total
    pct  = percent(frac)               # human‑readable percent
  )

# 2. Define your colors
my_cols <- c(
  "offs↑ only"  = "salmon",
  "diag↑ only"  = "plum3",
  "both↑"       = "lightsteelblue"
)

new_labels <- c(
  "offs↑ only" = "Predicted only by external location" ,
  "diag↑ only" = "Predicted only by single location" ,
  "both↑"      = "Predicted by both approaches"
)
### ---- links added to each pollinator by approach ----

thr <- best_discrete_threshold

diff_df_cat <- diff_df %>%
  mutate(
    # logical flags (robust to NA)
    diag_pred = !is.na(avg_sigm_predicted_diag) &
      avg_prop_isl_diag == 0 &
      avg_sigm_predicted_diag > thr,
    
    offs_pred = !is.na(avg_sigm_predicted_offs) &
      avg_prop_isl_diag == 0 &
      avg_sigm_predicted_offs > thr,
    
    pred_category = case_when(
      diag_pred & offs_pred ~ "both",
      diag_pred & !offs_pred ~ "local_only",
      !diag_pred & offs_pred ~ "external_only",
      TRUE ~ NA_character_
    )
  )

pollinator_link_counts <- diff_df_cat %>%
  filter(!is.na(pred_category)) %>%
  group_by(node_to, pred_category) %>%
  summarise(
    n_additional_links = n(),
    .groups = "drop"
  )

pollinator_link_counts_wide <- pollinator_link_counts %>%
  tidyr::pivot_wider(
    names_from  = pred_category,
    values_from = n_additional_links,
    values_fill = 0
  )

# sanity check
pollinator_link_counts_wide %>%
  mutate(total_predicted = local_only + external_only + both) %>%
  arrange(desc(total_predicted))

diff_df_cat %>%
  count(pred_category)

# show specific pollinators on map
highlight_poll <- "Camponotus_feae"

map_potential_missing_links_diag_offs <- ggplot(diff_df, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient) - we don't show existing links
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "white",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but predicted only with external data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_diag < best_discrete_threshold &
                                avg_sigm_predicted_offs > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_offs),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "salmon", high = "salmon",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # third layer: overlay only cells that were never observed but predicted only with within-system data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs < best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "plum", high = "plum",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # fourth layer: overlay interactions that were predicted using both approaches
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag == 0 & avg_sigm_predicted_offs > best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "lightsteelblue", high = "lightsteelblue",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_text(size = 4, angle = 90), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" ) +
  # ) + scale_x_discrete(
  #   labels = function(x) ifelse(x == "Camponotus_feae",
  #                               expression(italic("Camponotus feae")),
  #                               "")) + 
  tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  })) +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links_diag_offs)

pdf(
  file   = "map_missing_links_diags_offs_poll_names.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(map_potential_missing_links_diag_offs)
dev.off()     # close the file


### ---- Fig. 3b: pie chart ----
# 3. Make the pie
pie_chart <- ggplot(df_counts, aes(x = "", y = n, fill = sigm_cat)) +
  geom_col(width = 1, color = "white") +      # white border between slices
  coord_polar(theta = "y") +                  # convert bar → pie
  scale_fill_manual(values = my_cols,
                    labels = new_labels) +
  theme_void() +                              # remove axes/background
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position  = "bottom",
    legend.direction = "vertical"
  ) +
  #labs(title = "Interactions by Significance Category") +
  geom_text(
    aes(x = 1.2,label = n),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 6
  )

pie_chart
# pdf(
#   file   = "pie_chart.pdf",
#   width  = 7,    # inches
#   height = 7,
#   family = "Helvetica"   # or another installed font
# )
# print(pie_chart)
# dev.off()     # close the file
# 
# png(
#   filename = "pie_chart.png",
#   width    = 5,
#   height   = 5,
#   units    = "in",
#   res      = 300,
#   family   = "Helvetica"
# )
# 
# # --- your plotting code here ---
# print(pie_chart)
# 
# dev.off()
# 

### ---- only observed interactions ----
map_potential_missing_links_diag_offs_verified <- ggplot(diff_df, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient) - we don't show existing links
  geom_tile(aes(fill = avg_prop_isl_diag)) +
  scale_fill_gradient(low = "white", high = "white",
                      name = "Observed links:\nproportion\nof islands\nobserved",
                      breaks = seq(0, 1, 0.2)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but predicted only with external data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag != 0 & avg_sigm_predicted_diag < best_discrete_threshold &
                                avg_sigm_predicted_offs > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_offs),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "salmon", high = "salmon",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # third layer: overlay only cells that were never observed but predicted only with within-system data
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag != 0 & avg_sigm_predicted_offs < best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "plum", high = "plum",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # fourth layer: overlay interactions that were predicted using both approaches
  geom_tile(
    data = diff_df %>% filter(avg_prop_isl_diag != 0 & avg_sigm_predicted_offs > best_discrete_threshold &
                                avg_sigm_predicted_diag > best_discrete_threshold),
    aes(fill = avg_sigm_predicted_diag),
    alpha = 0.6
  ) +
  scale_fill_gradient(low = "lightsteelblue", high = "lightsteelblue",
                      name = "Predicted links:\naverage predicted\nprobability",
                      breaks = seq(0, 1, 0.1)) +
  
  # Final adjustments
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x = element_blank(), 
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.position = "bottom",         # Place legends at the bottom
    legend.box = "horizontal" 
  ) + tme +
  scale_y_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

print(map_potential_missing_links_diag_offs_verified)

map_potential_missing_links_diag_offs_verified_noleg <- map_potential_missing_links_diag_offs_verified +
  theme(legend.position = "none")

png(
  filename = "map_potential_missing_links_diag_offs_verified_noleg.png",
  width    = 12,
  height   = 6,
  units    = "in",
  res      = 300,
  family   = "Helvetica"
)

# --- your plotting code here ---
print(map_potential_missing_links_diag_offs_verified_noleg)

dev.off()

