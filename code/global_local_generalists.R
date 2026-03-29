# ---- correlation between global and local generalism ----

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


# mixed model
# Are species that are globally more generalized also more generalized locally, 
# accounting for repeated measures across islands?
library(lme4)

m1_poll <- lmer(island_degree ~ global_degree + 
             (1 | node_to),
           data = poll_degree_combined)

summary(m1_poll)

# variance in local generalization

poll_variation <- poll_degree_combined %>%
  group_by(node_to) %>%
  summarise(
    mean_local_degree = mean(island_degree),
    var_local_degree  = var(island_degree),
    sd_local_degree   = sd(island_degree),
    n_islands         = n(),
    .groups = "drop"
  )

hist(poll_variation$var_local_degree,
     breaks = 20,
     main = "Variability of local degree across islands (pollinators)",
     xlab = "Variation") + tme


plant_island_degree <- plant_degree_combined %>%
  mutate(node_from = fct_reorder(node_from, island_degree, .fun = median)) %>%
  ggplot(aes(x = node_from, y = island_degree, fill = node_from)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis") +
  # coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10,
                                   angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "Plant", y = "Island level degree") +
  tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

pdf(
  file   = "results/paper_figs/plant_island_degree.pdf",
  width  = 11,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
print(plant_island_degree)
dev.off()     # close the file

poll_island_degree <- poll_degree_combined %>%
  mutate(node_to = fct_reorder(node_to, island_degree, .fun = median)) %>%
  ggplot(aes(x = node_to, y = island_degree, fill = node_to)) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "magma") +
  # coord_flip() +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank()) + 
  labs(x = "Pollinator", y = "Island level degree") +
  tme +
  scale_x_discrete(labels = function(x) lapply(strsplit(x, "_"), function(y) {
    bquote(italic(.(paste(y, collapse = " "))))
  }))

pdf(
  file   = "results/paper_figs/poll_island_degree.pdf",
  width  = 11,    # inches
  height = 8,
  family = "Helvetica"   # or another installed font
)
print(poll_island_degree)
dev.off()     # close the file

# ---- correlate predicted and observed local degree ----
# if we want to see if the amount of added links (locally in the layer combination) is correlated with the local degree of the species:
df <- df %>%
  mutate(island_id = paste(train_layer, test_layer, sep = "_")) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values))

# plants
# average predicted probability across iterations
pred_avg <- df %>%
  group_by(train_layer, test_layer, node_from, node_to, removed) %>%
  summarise(
    mean_pred = mean(predicted_prob_sigm, na.rm = TRUE),
    original_links = first(original_links),  # constant across itr
    .groups = "drop"
  )

# predicted degree per species per layer combo
pred_degree <- pred_avg %>%
  filter(removed == 1,
         mean_pred > best_discrete_threshold) %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(
    predicted_degree = n_distinct(node_to),
    .groups = "drop"
  )

# observed local degree
obs_degree <- pred_avg %>%
  filter(original_links > 0) %>%
  group_by(train_layer, test_layer, node_from) %>%
  summarise(
    observed_degree = n_distinct(node_to),
    .groups = "drop"
  )

# combine
degree_comparison <- obs_degree %>%
  left_join(pred_degree,
            by = c("train_layer", "test_layer", "node_from")) %>%
  mutate(
    predicted_degree = replace_na(predicted_degree, 0)
  )

# correlate
correlation_plants_d <- cor.test(degree_comparison$observed_degree, degree_comparison$predicted_degree, use = "complete.obs", method = "pearson")
correlation_plants_d

# Extract correlation coefficient and p-value
r_value_d <- round(correlation_plants_d$estimate, 2)
p_value_d <- formatC(correlation_plants_d$p.value, digits = 2)  # or round as you prefer
label_text_plants_d <- paste0("r = ", r_value_d, ", p = ", p_value_d)

# plot 
plant_degree_obs_pred <- ggplot(degree_comparison, aes(x = observed_degree, y = predicted_degree)) +
  geom_point(alpha = 0.5, size = 2, color = "seagreen3") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Observed degree",
    y = "Predicted degree",
    title = paste("Plants:", label_text_plants_d)   # <--- add label in title
  ) +
  theme_minimal() + tme

plant_degree_obs_pred

# pollinators

# predicted degree per pollinator per layer combo (pollinators = node_to)
pred_degree_poll <- pred_avg %>%
  filter(removed == 1,
         mean_pred > best_discrete_threshold) %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(
    predicted_degree = n_distinct(node_from),
    .groups = "drop"
  )

# observed local degree per pollinator (FULL network; pollinators = node_to)
obs_degree_poll <- pred_avg %>%
  filter(original_links > 0) %>%
  group_by(train_layer, test_layer, node_to) %>%
  summarise(
    observed_degree = n_distinct(node_from),
    .groups = "drop"
  )

# combine
degree_comparison_poll <- obs_degree_poll %>%
  left_join(pred_degree_poll,
            by = c("train_layer", "test_layer", "node_to")) %>%
  mutate(
    predicted_degree = replace_na(predicted_degree, 0)
  )

# correlate
correlation_poll_d <- cor.test(
  degree_comparison_poll$observed_degree,
  degree_comparison_poll$predicted_degree,
  use = "complete.obs",
  method = "pearson"
)
correlation_poll_d

# Extract correlation coefficient and p-value
r_value_d <- round(correlation_poll_d$estimate, 2)
p_value_d <- formatC(correlation_poll_d$p.value, digits = 2)
label_text_poll_d <- paste0("r = ", r_value_d, ", p = ", p_value_d)

# plot
poll_degree_obs_pred <- ggplot(degree_comparison_poll, aes(x = observed_degree, y = predicted_degree)) +
  geom_point(alpha = 0.5, size = 2, color = "thistle") +
  geom_smooth(method = "lm", se = FALSE, color = "navy") +
  labs(
    x = "Observed degree",
    y = "Predicted degree",
    title = paste("Pollinators:", label_text_poll_d)
  ) +
  theme_minimal() + tme

poll_degree_obs_pred

# Remove individual axis labels
plant_degree_clean_d <- plant_degree_obs_pred +
  labs(x = NULL, y = NULL)

poll_degree_clean_d <- poll_degree_obs_pred +
  labs(x = NULL, y = NULL)

# Combine the two panels
main_panel <- plot_grid(
  plant_degree_clean_d,
  poll_degree_clean_d,
  ncol = 2,
  align = "hv"
)

# Add bottom (shared x-axis) label
with_x_label <- plot_grid(
  main_panel,
  ggdraw() + draw_label("Observed degree",
                        fontface = "bold",
                        size = 16),
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Add left (shared y-axis) label
final_plot_d <- plot_grid(
  ggdraw() + draw_label("Predicted degree (withheld links)",
                        angle = 90,
                        fontface = "bold",
                        size = 16),
  with_x_label,
  ncol = 2,
  rel_widths = c(0.08, 1)
)

final_plot_d

# # Save to PDF
pdf("results/paper_figs/local_degree_predicted_links.pdf", width = 10, height = 7)  # adjust size as needed
grid::grid.draw(final_plot_d)
dev.off()

# ---- performance by degree binning ----
# add degrees to prediction data frame
df2 <- df %>%
  left_join(obs_degree,
            by = c("train_layer", "test_layer", "node_from"))

# create degree binning
df2 <- df2 %>%
  mutate(degree_bin = ntile(observed_degree, 3))

# and label them
df2$degree_bin <- factor(df2$degree_bin,
                         levels = 1:3,
                         labels = c("specialists", "intermediate", "generalists"))

# evaluate
df3 <- df2 %>%
  mutate(predicted_binary = if_else(predicted_prob_sigm > best_discrete_threshold, 1, 0))

# for each iteration:
perf_iter <- df3 %>%
  filter(removed == 1) %>% 
  group_by(itr, degree_bin) %>%
  summarise(
    TP = sum(original_links > 0 & predicted_binary == 1),
    FN = sum(original_links > 0 & predicted_binary == 0),
    TN = sum(original_links == 0 & predicted_binary == 0),
    FP = sum(original_links == 0 & predicted_binary == 1),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    F05 = (1.25) * (precision * recall) / ((0.25 * precision) + recall),
    .groups = "drop"
  )

# average across iterations:
perf_summary <- perf_iter %>%
  group_by(degree_bin) %>%
  summarise(
    mean_F05 = mean(F05, na.rm = TRUE),
    sd_F05 = sd(F05, na.rm = TRUE)
  )

p <- kruskal.test(F05 ~ degree_bin, data = perf_iter)

degree_binning <- ggplot(perf_iter, aes(x = degree_bin, y = F05, fill = degree_bin)) +
  geom_boxplot(notch = TRUE) +
  scale_fill_brewer(palette = "Pastel2") + 
  labs(x = "Species degree class",
       y = expression(F[0.5]),
       title = paste0("Kruskal–Wallis p = ",
                      signif(p$p.value, 3)),
       y = expression(F[0.5])) +
  theme_minimal() +
  theme(legend.position = "none") + tme

pdf("results/paper_figs/degree_binning.pdf", width = 6, height = 6)  # adjust size as needed
grid::grid.draw(degree_binning)
dev.off()
