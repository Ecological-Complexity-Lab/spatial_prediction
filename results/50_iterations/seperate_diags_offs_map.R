# missing links map for diag and off-daigonal seperately

df_island_sep <- df_island %>%
  separate(island_id, into = c("island1", "island2"), sep = "_", convert = TRUE)

df_island_diag <- df_island_sep %>%
  filter(island1 == island2)

df_summary_diag <- df_island_diag %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

# 3. Remove interactions that were never observed and predicted as zero.
# df_summary <- df_summary %>% 
#   filter(!(avg_prop == 0 & avg_sigm_predicted < 0.5))

df_never_observed <- df_summary_diag %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_never_observed_poll <- df_summary_diag %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")

# for overall degree

never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")

never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")

# if we want to omit the Euphorbias
# never_plants_degree_overall <- never_plants_degree_overall %>%
#   filter(!node_from %in% c("Euphorbia_balsamifera_m", "Euphorbia_balsamifera_f"))

# correlation
df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 2)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, non-observed interactions",
    title = "Plants"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.2,   # Adjust depending on your data range
           label = label_text_plants,
           size = 5,
           color = "black")
plant_degree

df_to_correlate <- never_poll_degree_overall
df_to_correlate$x <- df_to_correlate$overall_poll_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_poll <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_poll
# Extract correlation coefficient and p-value
r_value <- round(correlation_poll$estimate, 2)
p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Overall degree",
    y = "Number of predicted, /nnon-observed interactions",
    title = "Pollinators"
  ) +
  theme_minimal() + tme +
  annotate("text",
           x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.2,   # Adjust depending on your data range
           label = label_text_polls,
           size = 5,
           color = "black")
poll_degree

final_plot <- combine_plots(plant_degree, poll_degree)

## ---- never-observed links ----
### ---- heatmap related to island proportion ----
# here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of islands in which these interactions were observed.

# order species by their degree
df_summary_diag <- df_summary_diag %>% left_join(overall_poll_degree, by="node_to")
df_summary_diag <- df_summary_diag %>% left_join(overall_plant_degree, by="node_from")

# Determine the order for plants based on overall_plant_degree
plant_order <- df_summary_diag %>%
  distinct(node_from, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(node_from)

# Determine the order for pollinators based on overall_poll_degree
poll_order <- df_summary_diag %>%
  distinct(node_to, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(node_to)

# Reset the levels for the species factors
df_summary_diag$node_from <- factor(df_summary_diag$node_from, levels = plant_order)
df_summary_diag$node_to   <- factor(df_summary_diag$node_to, levels = poll_order)


map_missing_links_diags <- ggplot(df_summary_diag, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop)) +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Proportion\nof islands\nobserved") +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = df_summary_diag %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
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

map_missing_links_diags

# now off-diagonals

df_island_offs <- df_island_sep %>%
  filter(island1 != island2)

df_summary_offs <- df_island_offs %>%
  group_by(node_from, node_to) %>%
  summarise(
    avg_prop = mean(observed, na.rm = TRUE),       # proportion of islands with observation
    avg_sigm_predicted = mean(island_sigm_predicted, na.rm = TRUE),
    n_islands = n(),  # number of islands contributing
    .groups = "drop"
  )

# 3. Remove interactions that were never observed and predicted as zero.
# df_summary <- df_summary %>% 
#   filter(!(avg_prop == 0 & avg_sigm_predicted < 0.5))

df_never_observed <- df_summary_offs %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_from) %>%
  summarise(count_never_observed = n(), .groups = "drop")

df_never_observed_poll <- df_summary_offs %>%
  filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold) %>%
  group_by(node_to) %>%
  summarise(count_never_observed = n(), .groups = "drop")

never_plants_degree <- df_never_observed %>% left_join(avg_plant_degree, by="node_from")

# for overall degree

never_plants_degree_overall <- df_never_observed %>% left_join(overall_plant_degree, by="node_from")

never_poll_degree_overall <- df_never_observed_poll %>% left_join(overall_poll_degree, by="node_to")

# if we want to omit the Euphorbias
# never_plants_degree_overall <- never_plants_degree_overall %>%
#   filter(!node_from %in% c("Euphorbia_balsamifera_m", "Euphorbia_balsamifera_f"))

# correlation
df_to_correlate <- never_plants_degree_overall
df_to_correlate$x <- df_to_correlate$overall_plant_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_plants <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_plants
# Extract correlation coefficient and p-value
r_value <- round(correlation_plants$estimate, 2)
p_value <- formatC(correlation_plants$p.value, digits = 2)  # or round as you prefer
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

# plant_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
#   geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
#   geom_smooth(method = "lm", se = FALSE, color = "navy") +
#   labs(
#     x = "Overall degree",
#     y = "Number of predicted, non-observed interactions",
#     title = "Plants"
#   ) +
#   theme_minimal() + tme +
#   annotate("text",
#            x = Inf,
#            y = Inf,
#            hjust = 1.1,
#            vjust = 1.2,   # Adjust depending on your data range
#            label = label_text_plants,
#            size = 5,
#            color = "black")
# plant_degree

df_to_correlate <- never_poll_degree_overall
df_to_correlate$x <- df_to_correlate$overall_poll_degree
df_to_correlate$y <- df_to_correlate$count_never_observed

correlation_poll <- cor.test(df_to_correlate$x, df_to_correlate$y, use = "complete.obs", method = "pearson")
correlation_poll
# Extract correlation coefficient and p-value
r_value <- round(correlation_poll$estimate, 2)
p_value <- formatC(correlation_poll$p.value, digits = 2)  # or round as you prefer
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

# poll_degree <- ggplot(df_to_correlate, aes(x = x, y = y)) +
#   geom_point(alpha = 0.6, size = 2, color = "thistle") +
#   geom_smooth(method = "lm", se = FALSE, color = "navy") +
#   labs(
#     x = "Overall degree",
#     y = "Number of predicted, /nnon-observed interactions",
#     title = "Pollinators"
#   ) +
#   theme_minimal() + tme +
#   annotate("text",
#            x = Inf,
#            y = Inf,
#            hjust = 1.1,
#            vjust = 1.2,   # Adjust depending on your data range
#            label = label_text_polls,
#            size = 5,
#            color = "black")
# poll_degree

# final_plot <- combine_plots(plant_degree, poll_degree)

## ---- never-observed links ----
### ---- heatmap related to island proportion ----
# here we visualize the links that were never observed yet predicted to exist by the algorithm, and alongside them interactions that were observed, and the proportion of islands in which these interactions were observed.

# order species by their degree
df_summary_offs <- df_summary_offs %>% left_join(overall_poll_degree, by="node_to")
df_summary_offs <- df_summary_offs %>% left_join(overall_plant_degree, by="node_from")

# Determine the order for plants based on overall_plant_degree
plant_order <- df_summary_offs %>%
  distinct(node_from, overall_plant_degree) %>%
  arrange(desc(overall_plant_degree)) %>%
  pull(node_from)

# Determine the order for pollinators based on overall_poll_degree
poll_order <- df_summary_offs %>%
  distinct(node_to, overall_poll_degree) %>%
  arrange(desc(overall_poll_degree)) %>%
  pull(node_to)

# Reset the levels for the species factors
df_summary_offs$node_from <- factor(df_summary_offs$node_from, levels = plant_order)
df_summary_offs$node_to   <- factor(df_summary_offs$node_to, levels = poll_order)


map_missing_links_offs <- ggplot(df_summary_offs, aes(x = node_to, y = node_from)) +
  # First layer: background heatmap for proportion observed (blue gradient)
  geom_tile(aes(fill = avg_prop)) +
  scale_fill_gradient(low = "white", high = "steelblue", 
                      name = "Proportion\nof islands\nobserved") +
  
  # Reset fill scale so the next layer can have its own gradient
  new_scale_fill() +
  
  # Second layer: overlay only cells that were never observed but have high predicted value
  geom_tile(
    data = df_summary_offs %>% filter(avg_prop == 0, avg_sigm_predicted > best_discrete_threshold),
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

map_missing_links_offs

# what's the difference?
identical(df_summary_offs$node_from, df_summary_diag$node_from)

# 1. Which links only appear in one of the two data-frames?
only_offs <- df_summary_offs %>%
  select(node_from, node_to) %>%
  anti_join(df_summary_diag, by = c("node_from","node_to")) %>%
  mutate(source = "offs")

only_diag <- df_summary_diag %>%
  select(node_from, node_to) %>%
  anti_join(df_summary_offs, by = c("node_from","node_to")) %>%
  mutate(source = "diag")

bind_rows(only_offs, only_diag) %>%
  arrange(node_from, node_to)
# → these are the “new” / “dropped” interactions

# 2. For the shared links, compute per-column differences.

diff_df <- full_join(
  df_summary_offs,
  df_summary_diag,
  by     = c("node_from","node_to"),
  suffix = c("_offs","_diag")
) %>%
  mutate(
    diff_avg_prop            = avg_prop_offs            - avg_prop_diag,
    diff_avg_sigm_predicted  = avg_sigm_predicted_offs  - avg_sigm_predicted_diag,
    diff_n_islands           = n_islands_offs           - n_islands_diag,
    diff_overall_poll_degree = overall_poll_degree_offs - overall_poll_degree_diag,
    diff_overall_plant_degree= overall_plant_degree_offs- overall_plant_degree_diag
  )

view(diff_df)

# these are interactions that were predicted to exist in the off-diagonals but not in the diagonal
diff_df_filtered <- diff_df %>% filter(avg_sigm_predicted_diag < 0.6 & avg_sigm_predicted_offs >= 0.6)

# and these are interactions that were predicted to exist in the diagonal but not in the off-diagonals
diff_df_diags <- diff_df %>% filter(avg_sigm_predicted_offs < 0.6 & avg_sigm_predicted_diag >= 0.6)

# # View only those with any non‐zero difference:
# diff_df %>%
#   filter(
#     diff_avg_prop           != 0 |
#       diff_avg_sigm_predicted != 0 |
#       diff_n_islands          != 0 |
#       diff_overall_poll_degree!= 0 |
#       diff_overall_plant_degree!=0
#   ) %>%
#   arrange(desc(abs(diff_avg_sigm_predicted)))

# plot the differences


df_plot <- diff_df %>%
  filter(avg_prop_diag == 0) %>% 
  mutate(
    sigm_cat = case_when(
      avg_sigm_predicted_diag < best_discrete_threshold &
        avg_sigm_predicted_offs >= best_discrete_threshold ~ "offs↑ only",
      
      avg_sigm_predicted_offs < best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "diag↑ only",
      
      avg_sigm_predicted_offs >= best_discrete_threshold &
        avg_sigm_predicted_diag >= best_discrete_threshold ~ "both↑",
      
      TRUE ~ NA_character_
    )
  )

# 2. build the plot
ggplot(df_plot, aes(x = node_to, y = node_from)) +
  
  # # background = observed proportion
  # geom_tile(aes(fill = avg_prop_offs)) +
  # scale_fill_gradient(
  #   low  = "white",
  #   high = "white",
  #   name = "Prop\nobserved"
  # ) +
  
  # allow a second fill scale
  new_scale_fill() +
  
  # overlay our three categories, with fixed colours
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
    name   = "Difference in \npredicted probability",
    labels = c(
      "offs↑ only" = "offs ≥ 0.6\n& diag < 0.6",
      "diag↑ only" = "diag ≥ 0.6\n& offs < 0.6",
      "both↑"       = "both ≥ 0.6"
    )
  ) +
  
  # tidy up
  theme_minimal() +
  labs(x = "Pollinator", y = "Plant") +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 8),
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

# how many links did each category add?
df_plot %>%
  filter(avg_prop_diag == 0) %>% 
  group_by(sigm_cat) %>%
  summarise(
    n_links = n()
  )
