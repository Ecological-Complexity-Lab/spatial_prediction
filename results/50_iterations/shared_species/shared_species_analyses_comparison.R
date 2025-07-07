# ---- comparing analyses quality ----
## ---- load libraries ----
library(dplyr)
library(ggplot2)
library(ggpubr)

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
result_summary_shared_plants <- read.csv('working_df_shared_plants_island_50_itr.csv')
result_summary_shared_pollinators <- read.csv('working_df_shared_pollinators_island_100_itr.csv')
  
# 1) bind your four tables, giving each an ID
df_all <- bind_rows(
  result_summary_all %>% mutate(dataset = "all"),
  result_summary_shared_species %>% mutate(dataset = "shared_species"),
  result_summary_shared_plants %>% mutate(dataset = "shared_plants"),
  result_summary_shared_pollinators %>% mutate(dataset = "shared_pollinators")
)

# 2) boxplot with overall ANOVA p-value
ggplot(df_all, aes(x = dataset, y = nnse)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  stat_compare_means(
    method    = "anova",           # overall ANOVA
    label     = "p.format",        # formatted p-value
    label.y   = max(df_all$nse)*1.05
  ) +
  labs(
    title = "Distribution of NSE across networks",
    x     = "Network subset",
    y     = "Normalized Nash–Sutcliffe Efficiency (nNSE)"
  ) +
  theme_minimal() + tme

# —— Optional: pairwise Wilcoxon tests —— #
comparisons <- list(
  c("all","shared_species"), c("all","shared_plants"), c("all","shared_pollinators"),
  c("shared_species","shared_plants"), c("shared_species","shared_pollinators"), c("shared_plants","shared_pollinators")
)
ggplot(df_all, aes(dataset, nnse)) +
  geom_boxplot(alpha = 0.7) +
  stat_compare_means(
    comparisons = comparisons,
    method      = "wilcox.test",
    label       = "p.signif"
  ) +
  labs(x="Network", y="nNSE") +
  theme_minimal() + tme
