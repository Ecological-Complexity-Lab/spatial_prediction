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


tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())

## ---- plot ----
result_summary_isl <- read.csv('working_df_island_distance_fidelity_jaccard_netsize.csv')

layer_to_layer_plot_all_itr_isl <- 
  ggplot(result_summary_isl, aes(x = train_layer_name, y = test_layer_name, fill = f1_score)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary_isl[result_summary_isl$train_layer == result_summary_isl$test_layer, ],
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

print(layer_to_layer_plot_all_itr_isl)

result_summary_site <- read.csv('working_df_all_itr_60_binary_names.csv')

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

print(layer_to_layer_plot_all_itr_site)

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


