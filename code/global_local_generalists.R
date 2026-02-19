# correlation between global and local generalism

df_itr1_observed <- df %>% filter(itr == 1 & original_links > 0)

# plants
plant_island_degree <- df_itr1_observed %>%
  group_by(test_layer, node_from) %>%
  summarise(island_degree = n_distinct(node_to), .groups = "drop")

plant_global_degree <- df_itr1_observed %>%
  group_by(node_from) %>%
  summarise(global_degree = n_distinct(node_to), .groups = "drop")

plant_degree_combined <- plant_island_degree %>%
  left_join(plant_global_degree, by = "node_from")

cor.test(
  plant_degree_combined$island_degree,
  plant_degree_combined$global_degree,
  method = "pearson"   
)

# pollinators
poll_island_degree <- df_itr1_observed %>%
  group_by(test_layer, node_to) %>%
  summarise(island_degree = n_distinct(node_from), .groups = "drop")

poll_global_degree <- df_itr1_observed %>%
  group_by(node_to) %>%
  summarise(global_degree = n_distinct(node_from), .groups = "drop")

poll_degree_combined <- poll_island_degree %>%
  left_join(poll_global_degree, by = "node_to")

cor.test(
  poll_degree_combined$island_degree,
  poll_degree_combined$global_degree,
  method = "pearson"   
)

# plot
# plants
df_to_correlate_plants <- plant_degree_combined %>%
  mutate(
    x = global_degree,
    y = island_degree
  )

correlation_plants <- cor.test(
  df_to_correlate_plants$x,
  df_to_correlate_plants$y,
  use = "complete.obs",
  method = "pearson" # or "spearman"
)

r_value <- round(unname(correlation_plants$estimate), 2)
p_value <- formatC(correlation_plants$p.value, digits = 2, format = "g")
label_text_plants <- paste0("r = ", r_value, ", p = ", p_value)

plant_degree <- ggplot(df_to_correlate_plants, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Global degree",
    y = "Island-level degree",
    title = paste("Plants:", label_text_plants)
  ) +
  theme_minimal() + tme

# pollinators
df_to_correlate_poll <- poll_degree_combined %>%
  mutate(
    x = global_degree,
    y = island_degree
  )

correlation_poll <- cor.test(
  df_to_correlate_poll$x,
  df_to_correlate_poll$y,
  use = "complete.obs",
  method = "pearson" # or "spearman"
)

r_value <- round(unname(correlation_poll$estimate), 2)
p_value <- formatC(correlation_poll$p.value, digits = 2, format = "g")
label_text_polls <- paste0("r = ", r_value, ", p = ", p_value)

poll_degree <- ggplot(df_to_correlate_poll, aes(x = x, y = y)) +
  geom_point(alpha = 0.6, size = 2, color = "rosybrown2") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Global degree",
    y = "Island-level degree",
    title = paste("Pollinators:", label_text_polls)
  ) +
  theme_minimal() + tme

# combine
plant_degree_clean <- plant_degree + labs(x = NULL, y = NULL)
poll_degree_clean  <- poll_degree  + labs(x = NULL, y = NULL)

final_plot <- plot_grid(
  plot_grid(plant_degree_clean, poll_degree_clean, ncol = 2, align = "hv"),
  ggdraw() + draw_label("Global degree", fontface = "bold", size = 16),
  ncol = 1,
  rel_heights = c(1, 0.08)
)

final_plot <- plot_grid(
  ggdraw() + draw_label("Island-level degree",
                        angle = 90, fontface = "bold", size = 16),
  final_plot,
  ncol = 2,
  rel_widths = c(0.08, 1)
)

final_plot

pdf(
  file   = "results/paper_figs/global_local_degree.pdf",
  width  = 8,    # inches
  height = 5,
  family = "Helvetica"   # or another installed font
)
print(final_plot)
dev.off()     # close the file

