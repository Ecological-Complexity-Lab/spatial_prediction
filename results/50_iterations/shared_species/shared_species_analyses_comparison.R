# ---- comparing analyses quality ----
## ---- load libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(VennDiagram)


## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- load data ----
result_summary_all <- read.csv('working_df_all_species_island_50_itr.csv')
result_summary_shared_species <- read.csv('working_df_shared_species_island_50_itr.csv')
result_summary_shared_plants <- read.csv('working_df_shared_plants_island_100_itr.csv')
result_summary_shared_pollinators <- read.csv('working_df_shared_pollinators_island_100_itr.csv')
  
# 1) bind your four tables, giving each an ID
df_all <- bind_rows(
  result_summary_all %>% mutate(dataset = "all"),
  result_summary_shared_species %>% mutate(dataset = "shared_species"),
  result_summary_shared_plants %>% mutate(dataset = "shared_plants"),
  result_summary_shared_pollinators %>% mutate(dataset = "shared_pollinators")
)

## ---- plot distances between analyses ----
# 2) boxplot with overall ANOVA p-value
ggplot(df_all, aes(x = dataset, y = nnse)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(
    method    = "anova",           # overall ANOVA
    label     = "p.format",        # formatted p-value
    label.y   = max(df_all$nnse)*1.05
  ) +
  labs(
    title = "Distribution of NSE across networks",
    x     = "Network subset",
    y     = "Normalized Nash–Sutcliffe Efficiency (nNSE)"
  ) +
  theme_minimal() + tme

ggplot(df_all, aes(x = dataset, y = f1_score)) +
  geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(
    method    = "anova",           # overall ANOVA
    label     = "p.format",        # formatted p-value
    label.y   = max(df_all$f1_score)*1.05
  ) +
  labs(
    #title = "Distribution of NSE across networks",
    x     = "Network subset",
    y     = "F1 score"
  ) +
  theme_minimal() + tme

# —— Optional: pairwise Wilcoxon tests —— #
comparisons <- list(
  c("all","shared_species"), c("all","shared_plants"), c("all","shared_pollinators"),
  c("shared_species","shared_plants"), c("shared_species","shared_pollinators"), c("shared_plants","shared_pollinators")
)

ggplot(df_all, aes(dataset, nnse)) +
  geom_boxplot(notch = TRUE, alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.signif"
  ) +
  #scale_fill_manual(values = c("all" = "lightsteelblue", "shared_species" = "orange", "shared_plants" = "lightseagreen", "shared_pollinators" = "thistle")) +  # Assign custom colors to groups
  labs(x="Network subset", y="nNSE") +
  theme_minimal() + tme

## ---- venn diagram ----
# load more data...
missing_links_all_species <- read.csv("predicted_non_observed_links.csv", row.names = NULL)
missing_links_shared_species <- read.csv("predicted_non_observed_links_shared_species.csv", row.names = NULL)
missing_links_shared_plants <- read.csv("predicted_non_observed_links_plants.csv", row.names = NULL)
missing_links_shared_pollinators <- read.csv("predicted_non_observed_links_shared_pollinators.csv", row.names = NULL)

# Example: Suppose you have four data frames: df1, df2, df3, df4
# You will combine `node_to` and `node_from` to create unique interaction identifiers for each data frame

# Combine 'node_to' and 'node_from' columns for each data frame to create interaction sets
all_species_interactions <- unique(paste(missing_links_all_species$node_to, missing_links_all_species$node_from, sep = "_"))
shared_species_interactions <- unique(paste(missing_links_shared_species$node_to, missing_links_shared_species$node_from, sep = "_"))
shared_plants_interactions <- unique(paste(missing_links_shared_plants$node_to, missing_links_shared_plants$node_from, sep = "_"))
shared_pollinators_interactions <- unique(paste(missing_links_shared_pollinators$node_to, missing_links_shared_pollinators$node_from, sep = "_"))

# Create a list of these sets to pass to the Venn diagram function
interaction_sets <- list(
  "All species" = all_species_interactions,
  "Shared species" = shared_species_interactions,
  "Shared plants" = shared_plants_interactions,
  "Shared pollinators" = shared_pollinators_interactions
)

# Plot the Venn diagram
venn.plot <- venn.diagram(
  x = interaction_sets,
  category.names = c("All species", "Shared species", "Shared plants", "Shared pollinators"),
  filename = NULL,  # You can save it as a file if needed
  output = TRUE,
  col = "transparent",  # Outline color of the circles
  fill = c("lightsteelblue", "orange", "lightseagreen", "thistle"),  # Colors for each circle
  alpha = 0.5,  # Transparency level for fill color
  cex = 1.5,  # Text size
  fontface = "bold",  # Text boldness
  fontfamily = "sans",  # Text font
  cat.cex = 1.5,  # Text size for the category names
  cat.fontface = "bold",  # Text font for category names
  cat.pos = 0,  # Positioning of the category names (0 means below)
  cat.dist = 0.1  # Distance of category names from circles
)

# Clear the plot window before drawing the Venn diagram
grid.newpage()  # Clear any existing plot

# Draw the Venn diagram
grid.draw(venn.plot)
