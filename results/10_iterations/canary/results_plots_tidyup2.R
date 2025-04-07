## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)
library(emln)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)  # for get_legend()
library(patchwork)


tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- heatmaps ----
result_summary_isl <- read.csv('working_df_island_distance_fidelity_jaccard_netsize.csv')
result_summary_isl <- read.csv('result_summary_canary_scaled_with_distance.csv') # scaled
result_summary_site <- read.csv('working_df_all_itr_60_binary_scaled_site_names.csv') # scaled, site
result_summary_site <- read.csv('working_df_all_itr_60_weighted_scaled_site.csv') # weighted scaled site

layer_to_layer_plot_all_itr_isl <- 
  ggplot(result_summary_site, aes(x = train_layer_name, y = test_layer_name, fill = recall)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "black", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_site[result_summary_site$train_layer == result_summary_site$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +  # Set NA values to gray
  labs(x = "Added layer", y = "Predicted layer", fill = "Recall") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed() + tme

print(layer_to_layer_plot_all_itr_isl)

result_summary_site <- read.csv('working_df_all_itr_60_weighted_scaled_site_names.csv')

layer_to_layer_plot_all_itr_site <- 
  ggplot(result_summary_site, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_site[result_summary_site$train_layer == result_summary_site$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "f1 score") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

layer_to_layer_plot_all_itr_site_pr <- 
  ggplot(result_summary_site, aes(x = train_layer_name, y = test_layer_name, fill = balanced_accuracy)) +
  geom_tile(color = "black", linewidth = 0.1) +  
  geom_tile(data = result_summary_site[result_summary_site$train_layer == result_summary_site$test_layer, ],
            color = "black", linewidth = 1.2) +
  scale_fill_gradient2(low = "steelblue2", mid = "white", high = "salmon2", 
                       midpoint = 0.5, na.value = "gray") +
  labs(x = "Added layer", y = "Predicted layer", fill = "Balanced \naccuracy") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  coord_fixed() + tme


print(layer_to_layer_plot_all_itr_site_pr)

# --- Remove individual legends and axis titles from the plots ---
# Reduce margins to decrease extra spacing between heatmaps.
p1 <- layer_to_layer_plot_all_itr_site +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

p2 <- layer_to_layer_plot_all_itr_isl +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

# --- Arrange the two heatmaps with adjusted spacing ---
# Here we set equal relative widths for p1 and p2.
combined_plots <- arrangeGrob(
  p1, p2, 
  ncol = 2, 
  widths = c(1, 1)
)

# --- Add common overall axis labels ---
# Adjust the x-axis label’s vertical position by changing vjust.
combined_with_axes <- arrangeGrob(
  combined_plots,
  bottom = textGrob("Training layer", gp = gpar(fontsize = 14), vjust = -1.5),
  left   = textGrob("Predicted layer", rot = 90, gp = gpar(fontsize = 14))
)

# --- Extract a common legend with adjusted key height ---
legend_plot <- layer_to_layer_plot_all_itr_site +
  theme(legend.position = "right",
        legend.key.height = unit(0.5, "cm"))
common_legend <- get_legend(legend_plot)

# --- Center the legend vertically ---
# Wrapping the legend in arrangeGrob with heights set to unit(1, "npc") makes it fill its column height.
legend_centered <- arrangeGrob(common_legend, heights = unit(1, "npc"))

# --- Arrange the combined heatmaps and legend side by side ---
# Adjust widths to control the relative space; here the legend column is narrower.
final_plot <- grid.arrange(
  combined_with_axes,
  legend_centered,
  ncol = 2,
  widths = c(2, 0.3)
)

print(final_plot)

# multivariable analysis
library(GGally)
# Select variables of interest including f1_score and your explanatory variables
selected_vars <- result_summary_isl[, c("f1_score", "distance_km", "avg_sorensen_plants",	"avg_sorensen_pollinators",	"jaccard_pollinators",	"jaccard_plants",	"jaccard_edges",	"size_P",	"density_P",	"size_C",	"density_C"
)]

# Create a pairwise scatterplot matrix
ggpairs(selected_vars) + tme

# ---- distance ----
### ---- calculate distance between islands ----
# first convert distances to distances between islands
island_scale <- read.csv('working_df_evaluators_island_names.csv')
island_scale <- read.csv('working_df_all_itr_60_binary_scaled_island.csv') # scaled version
island_scale <- read.csv('working_df_all_itr_60_weighted_scaled_island.csv') # weighted, scaled version

distance_table <- read.csv("distance_between_sites_canary.csv", row.names = NULL)

# Function to extract island names (removes "_site_X")
extract_island <- function(name) {
  gsub("_site_[12]", "", name)
}

# Create new table with averaged distances at the island level
distance_island_table <- distance_table %>%
  mutate(
    from_island = extract_island(from),
    to_island = extract_island(to)
  ) %>%
  group_by(from_island, to_island) %>%
  summarise(
    avg_distance_m = mean(distance_m),
    avg_distance_km = mean(distance_km),
    .groups = "drop"
  ) %>%
  mutate(
    avg_distance_m = ifelse(from_island == to_island, 0, avg_distance_m),
    avg_distance_km = ifelse(from_island == to_island, 0, avg_distance_km)
  ) %>%
  rename(from = from_island, to = to_island)  # Rename after calculation

# Print result
print(distance_island_table)

# Modify the 'from' and 'to' columns in distance_island_table
distance_island_table <- distance_island_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

# Add names and distances to the main table
net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Step 1: Create a new grouped tibble
new_layer_names <- net_name %>%
  mutate(group_id = (layer_id + 1) %/% 2) %>%  # Group pairs into 1, 2, 3...
  group_by(group_id) %>%
  summarise(name = gsub(" site.*", "", first(name)), .groups = "drop")  # Keep only location name

# Add to main table
island_scale <- island_scale %>%
  left_join(new_layer_names, by = c("train_layer" = "group_id")) %>%
  rename(train_layer_name = name) %>%
  left_join(new_layer_names, by = c("test_layer" = "group_id")) %>%
  rename(test_layer_name = name)

# Add distances
island_scale <- island_scale %>%
  left_join(distance_island_table, by = c("train_layer_name" = "from", "test_layer_name" = "to")) %>%
  mutate(distance_km = avg_distance_km)

island_scale <- island_scale %>% select(-avg_distance_km)

write.csv(island_scale, 'working_df_island_weighted_scaled_evaluators_distance.csv')

# for site scale
site_scale <- read.csv('working_df_all_itr_60_binary_scaled_site_names.csv') # scaled version

# Modify the 'from' and 'to' columns in distance_island_table
distance_table <- distance_table %>%
  mutate(from = gsub("_", " ", from),
         to = gsub("_", " ", to))

# Add to main table

site_scale <- site_scale %>%
  left_join(
    distance_table,
    by = c("train_layer_name" = "from", "test_layer_name" = "to")
  ) %>%
  mutate(distance_km = if_else(train_layer_name == test_layer_name,
                               0,              # distance = 0 if same site
                               distance_km))   # otherwise, keep joined distance

write.csv(site_scale, 'working_df_site_scaled_evaluators_distance.csv')

### ---- correlate with distance ----
distance_site <- read.csv('result_summary_canary_with_distance.csv')
distance_site <- read.csv('working_df_site_scaled_evaluators_distance.csv') # scaled, binary
result_summary_isl <- read.csv('working_df_islands_scaled_evaluators_distance.csv') # scaled, binary
distance_site <- read.csv('result_netsize_canaries_distance_names_site_weighted_scaled.csv') # scaled, weighted
result_summary_isl <- read.csv('working_df_island_weighted_scaled_evaluators_distance.csv') # scaled, weighted

correlation_site <- cor.test(distance_site$recall, distance_site$distance_km, use = "complete.obs", method = "pearson")
correlation_site
# Extract correlation coefficient and p-value
r_value_site <- round(correlation_site$estimate, 3)
p_value_site <- formatC(correlation_site$p.value, format = "f", digits = 2)

correlation_island <- cor.test(result_summary_isl$recall, result_summary_isl$distance_km, use = "complete.obs", method = "pearson")
correlation_island
# Extract correlation coefficient and p-value
r_value_island <- round(correlation_island$estimate, 3)
p_value_island <- formatC(correlation_island$p.value, format = "f", digits = 2)

label_text_site <- paste0("r = ", r_value_site, ", p = ", p_value_site)
label_text_isl <- paste0("r = ", r_value_island, ", p = ", p_value_island)

cor_plot_site <- 
  ggplot(distance_site, aes(x = distance_km, y = recall)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "Recall") +
  annotate("text",
           x = 365, y = 1.05,   # Adjust depending on your data range
           label = label_text_site,
           size = 3,
           color = "black") +
  tme

cor_plot_isl <- 
  ggplot(result_summary_isl, aes(x = distance_km, y = recall)) +
  geom_point(color = "salmon2", size = 2) +  # Points representing pairs of layers
  geom_smooth(method = "lm", se = FALSE, color = "steelblue2") +  # Linear regression line
  labs(x = "Geographical distance (km)",
       y = "Recall") +
  annotate("text",
           x = 362, y = 1.05,   # Adjust depending on your data range
           label = label_text_isl,
           size = 3,
           color = "black") +
  tme

p1 <- cor_plot_site +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 0.3), "cm"))

p2 <- cor_plot_isl +
  theme(legend.position = "none",
        axis.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 1, 0.3), "cm"))

combined_plots <- arrangeGrob(
  p1, p2, 
  ncol = 2, 
  widths = c(1, 1)
)
combined_with_axes <- arrangeGrob(
  combined_plots,
  bottom = textGrob("Geographical distance (km)", gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
  left   = textGrob("Recall", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

final_plot <- grid.arrange(
  combined_with_axes,
  ncol = 2,
  widths = c(2, 0.3)
)

## ---- distribution ----

isl_distrib <- ggplot(result_summary_isl, aes(x = f1_score)) +
  geom_histogram(bins = 20, fill = "rosybrown2", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(title = "Island scale",
       x = "F1 score",
       y = "Count") +
  tme

site_specificity <- ggplot(result_summary_site, aes(x = specificity)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Specificity",
       y = "Count") +
  tme

site_f1 <- ggplot(result_summary_site, aes(x = f1_score)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "F1 score",
       y = "Count") +
  tme

site_ba <- ggplot(result_summary_site, aes(x = balanced_accuracy)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Balanced accuracy",
       y = "Count") +
  tme

site_precision <- ggplot(result_summary_site, aes(x = precision)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Precision",
       y = "Count") +
  tme

site_recall <- ggplot(result_summary_site, aes(x = recall)) +
  geom_histogram(bins = 20, fill = "lightsteelblue", color = "black", alpha = 0.5) + 
  #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
  labs(x = "Recall",
       y = "Count") +
  tme

p1 <- site_f1 +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p2 <- site_recall +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p3 <- site_ba +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p4 <- site_precision +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

p5 <- site_specificity +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.1, 0.3), "cm"))

combined_plots <- arrangeGrob(
  p1, p2, p3, p4, p5,
  ncol = 3, 
  nrow = 2
)
combined_with_axes <- arrangeGrob(
  combined_plots,
  #bottom = textGrob("F1 score", gp = gpar(fontsize = 14, fontface = "bold"), vjust = -1.5),
  left   = textGrob("Count of instances", rot = 90, gp = gpar(fontsize = 14, fontface = "bold"))
)

final_plot <- grid.arrange(
  combined_with_axes,
  ncol = 3,
  widths = c(2, 0.3, 0.3)
)

# ---- diagonal and offs comparison ----
result_table_site <- read.csv('working_df_all_itr_60_binary_names.csv') # site scale
result_table_site <- read.csv('working_df_all_itr_60_binary_scaled_site_names.csv') #scaled, site

result_table_site <- result_table_site %>%
  mutate(layer_comparison = case_when(
    train_layer == test_layer ~ "Diagonal",
    train_layer != test_layer ~ "Off-diagonals"
  ))

custom_colors <- c("Diagonal" = "steelblue",
                   "Off-diagonals" = "thistle")

plot_boxplot <- function(data, metric, y_axis_label = "Balanced accuracy", 
                            stat_label_y = NULL, stat_size = 3) {
  ggplot(data, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = FALSE, alpha = 0.4, color = "black") +
    theme_minimal() +
    labs(y = y_axis_label) +   # y-axis title is set via the function argument
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    tme + 
    scale_fill_manual(values = custom_colors) +
    stat_compare_means(method = "t.test", label = "p.signif", hide.ns = FALSE, 
                       comparisons = list(c("Diagonal", "Off-diagonals")),
                       label.y = stat_label_y,  # Adjust vertical position here
                       size = stat_size) +
    scale_y_continuous(limits = c (0.3, 0.9))
    
}

plot_hist <- function(data, metric, 
                                        x_axis_label = "Balanced accuracy", 
                                        y_axis_label = "Count") {
  ggplot(data, aes(x = .data[[metric]], fill = layer_comparison)) +
    geom_histogram(aes(y = ..count..), alpha = 0.4, color = "black", bins = 8, position = "dodge") +
    #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
    theme_minimal() +
    labs(x = x_axis_label,
         y = y_axis_label,
         fill = "Layer comparison") +
    theme(
      axis.text.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(values = custom_colors) + tme
}

# if you prefer density over counts:
plot_hist_density <- function(data, metric, 
                      x_axis_label = "Balanced accuracy", 
                      y_axis_label = "Density") {  
  ggplot(data, aes(x = .data[[metric]], fill = layer_comparison)) +
    geom_histogram(aes(y = after_stat(density)), alpha = 0.4, color = "black", bins = 8, position = "dodge") +
    #geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", linewidth = 1) +
    theme_minimal() +
    labs(x = x_axis_label,
         y = y_axis_label,
         fill = "Layer comparison") +
    theme(
      axis.text.x = element_text(hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    scale_fill_manual(values = custom_colors) + tme
}

# boxplots:
site_ba <- plot_boxplot(result_table_site, metric = "balanced_accuracy", 
                           y_axis_label = "Balanced accuracy", stat_label_y = 0.85, stat_size = 6)
site_f1 <- plot_boxplot(result_table_site, metric = "f1_score", 
                           y_axis_label = "F1 score", stat_label_y = 0.85, stat_size = 6)

combined_plots <- arrangeGrob(
  site_ba, site_f1, 
  ncol = 2, 
  widths = c(1, 1)
)

final_plot <- grid.arrange(
  combined_plots,
  ncol = 2,
  widths = c(2, 0.3)
)

# histograms: 
hist_ba <- plot_hist(result_table_site, metric = "balanced_accuracy", 
                        y_axis_label = "Count")
hist_f1 <- plot_hist(result_table_site, metric = "f1_score", 
                        y_axis_label = "Count",
                        x_axis_label = "F1 score")

hist_f1 <- hist_f1 + theme(axis.title.y = element_blank())

combined_plot <- hist_ba + hist_f1 + 
  plot_layout(guides = "collect") +
  # Optionally, set the legend position (e.g., to the right or bottom)
  plot_annotation(theme = theme(legend.position = "right"))
# hist_ba <- hist_ba +
#   theme(legend.position = "none",
#         axis.title = element_blank())
# 
# combined_plots <- arrangeGrob(
#   hist_ba, hist_f1, 
#   ncol = 2, 
#   widths = c(1, 1)
# )
# 
# final_plot <- grid.arrange(
#   combined_plots,
#   ncol = 2,
#   widths = c(7, 0.4)
# )
