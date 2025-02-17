# ---- comparing scales ----
## ---- load libraries ----
library(dplyr)
library(tidyverse)
## ---- load data ----
site_scale <- read_csv('result_summary_canary_with_distance.csv') # site scale
island_scale # island scale

## ---- summarize predictive performance ----
# summarize max, min, average, se of f1, ba, recall for both scales

#plot
# Define custom colors for each group
custom_colors <- c("Diagonal" = "#1b9e77", 
                   "Off-diagonal (train > test)" = "#d95f02", 
                   "Off-diagonal (train < test)" = "#7570b3")

# Function to create notched boxplots with customizations
plot_f1_boxplot <- function(metric, y_axis_label = "F1 score") {
  ggplot(result_table, aes(x = layer_comparison, y = .data[[metric]], fill = layer_comparison)) +
    geom_boxplot(notch = TRUE, alpha = 0.4, color = "black") +  # Notched, semi-transparent, black outline
    theme_minimal() +
    labs(x = "Layer comparison",
         y = y_axis_label) +  # Custom y-axis title
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # Black frame around the plot
    ) +
    tme +
    scale_fill_manual(values = custom_colors) +  # Apply custom colors
    stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(
      c("Diagonal", "Off-diagonal (train > test)"),
      c("Diagonal", "Off-diagonal (train < test)"),
      c("Off-diagonal (train > test)", "Off-diagonal (train < test)")
    ))  # Pairwise t-tests with significance labels
}

# Generate boxplots for each F1 score metric
plot_f1_boxplot("f1_score")
