# ---- exploring results for the weighted predictions ----
## ---- themes ----
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 12, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- load libraries ----

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}
## ---- load data ----

nonbinary <- read_csv('nonbinary_equal_0_1_removal_60_0.csv')

nonbinary <- nonbinary %>% mutate(sigm_predicted = sigmoid(predicted_values))
nonbinary_filtered <- nonbinary %>% filter(k == 2) %>% filter(lambda == 0.1)
write.csv(nonbinary_filtered, 'working_all_itr_nonbinary_canary.csv')

nonbinary_removed <- nonbinary_filtered %>% filter(removed == 1)
# check relationship between predicted and observed values

correlation <- cor.test(nonbinary_removed$sigm_predicted, nonbinary_removed$original_links, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

ggplot(nonbinary_removed, aes(x = original_links, y = sigm_predicted)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  labs(x = "Weight of original links", y = "Predicted value (logistic transformation)") +
  theme_minimal() +
  tme + 
  annotate("text",
           x = 75, y = 0.9,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black") +
  theme(
    # Add a black frame around facet labels with thickness
    #strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black")
  )

## ---- correlation ----
# Calculate correlations per train_layer, test_layer, and iteration
correlation_results <- nonbinary_removed %>%
  group_by(train_layer, test_layer, itr) %>%
  summarise(correlation = cor(original_links, sigm_predicted, use = "complete.obs"), .groups = "drop")

# Compute average correlation per layer combination
average_correlation <- correlation_results %>%
  group_by(train_layer, test_layer) %>%
  summarise(avg_correlation = mean(correlation, na.rm = TRUE), .groups = "drop")

layer_to_layer_plot_corr <- 
  ggplot(average_correlation, aes(x = factor(train_layer, levels = unique(train_layer)),  
                                      y = factor(test_layer, levels = unique(test_layer)),
                                                 fill = avg_correlation)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = average_correlation[average_correlation$train_layer == average_correlation$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "correlation") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

## ---- distance between predicted and observed values ----
# Compute relative weight
nonbinary_filtered <- nonbinary_filtered %>%
  group_by(train_layer, test_layer, itr) %>%
  mutate(total_weight = sum(original_links, na.rm = TRUE),  # Total weight per group
         relative_weight = original_links / total_weight) %>%  # Relative weight calculation
  ungroup()

# View the updated dataframe
head(nonbinary_filtered)

nonbinary_filtered <- nonbinary_filtered %>% mutate(pred_obs = sigm_predicted - relative_weight)

## ---- add partner fidelity ----
# filter out cases in which train = test layer, and only occurring interactions
nonbinary_fidelity <- nonbinary_filtered %>% filter(train_layer < test_layer | train_layer == test_layer) %>% 
  filter(original_links == 1)

## ---- plants ----
# 1. For each plant (node_from) and layer, gather the pollinators (node_to).
#    (Assuming 'train_layer' is the relevant layer ID—adapt as needed.)
df_plant_partners <- nonbinary_fidelity %>%
  distinct(node_from, train_layer, node_to) %>%
  group_by(node_from, train_layer) %>%
  summarise(partners = list(unique(node_to)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen <- df_plant_partners %>%
  group_by(node_from) %>%
  summarise(
    mean_sorensen_plants = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

nonbinary_fidelity_merged <- nonbinary_fidelity %>%
  left_join(df_sorensen, by = "node_from")

# for pollinators

df_pollinators <- nonbinary_fidelity %>%
  distinct(node_to, train_layer, node_from) %>%
  group_by(node_to, train_layer) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# 2. For each plant, compute mean Sorensen similarity across all pairs of layers.
df_sorensen_pollinators <- df_pollinators %>%
  group_by(node_to) %>%
  summarise(
    mean_sorensen_pollinators = {
      n_layers <- n()
      # If the plant is only in one layer, there are no pairs, so return NA
      if (n_layers < 2) {
        NA_real_
      } else {
        # Get all pairwise combinations of rows in this group
        idx_pairs <- combn(n_layers, 2)
        # Compute Sorensen for each pair
        sims <- apply(idx_pairs, 2, function(idx) {
          p1 <- partners[[idx[1]]]
          p2 <- partners[[idx[2]]]
          a  <- length(intersect(p1, p2))          # shared partners
          b  <- length(setdiff(p1, p2))            # unique to first layer
          c  <- length(setdiff(p2, p1))            # unique to second layer
          2 * a / (2 * a + b + c)
        })
        mean(sims)  # average Sorensen for that plant
      }
    }
  )

nonbinary_fidelity_merged <- nonbinary_fidelity_merged %>%
  left_join(df_sorensen_pollinators, by = "node_to")

# write.csv(nonbinary_fidelity_merged, 'nonbinary_canaries_fidelity_obs_pred.csv')
nonbinary_fidelity_merged <- read.csv('nonbinary_canaries_fidelity_obs_pred.csv')

# plot relationship with partner fidelity
nonbinary_fidelity_merged_removed <- nonbinary_fidelity_merged %>% filter(removed == 1)
#nonbinary_fidelity_merged_removed_itr_1 <- nonbinary_fidelity_merged %>% filter(itr == 1) # try 1 itr
correlation <- cor.test(nonbinary_fidelity_merged_removed$pred_obs, nonbinary_fidelity_merged_removed$mean_sorensen_plants, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

ggplot(nonbinary_fidelity_merged_removed, aes(x = mean_sorensen_plants, y = pred_obs)) +
  geom_point(color = "seagreen3", alpha = 0.6, size = 2) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +  # Trendline
  labs(x = "Mean Sorensen similarity",
       y = "Predicted - observed",
       title = "Sigmoid predicted value - relative observed weight vs. mean partner fidelity of plants (Sorensen)") +
  tme +
  annotate("text",
           x = 0.92, y = 1.05,   # Adjust depending on your data range
           label = label_text,
           size = 3.5,
           color = "black")

# plotting them together

# Reshape data: Convert wide to long format
long_data <- nonbinary_fidelity_merged_removed %>%
  pivot_longer(cols = c(mean_sorensen_pollinators, mean_sorensen_plants),
               names_to = "group", values_to = "sorensen_value")

# Plot both on the same graph with different colors
ggplot(long_data, aes(x = sorensen_value, y = pred_obs, color = group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("mean_sorensen_pollinators" = "thistle", "mean_sorensen_plants" = "aquamarine4"),
                     labels = c("Pollinators", "Plants")) +
  labs(x = "Mean Sorensen similarity",
       y = "Predicted - observed",
       title = "Predicted - Observed vs. Mean Partner Fidelity",
       color = "Group") +
  theme_minimal() +
  tme

# plot the distribution of partner fidelity with a histogram
nonbinary_fidelity_merged_removed_itr_1_tn_1_ts_2 <- nonbinary_fidelity_merged_removed_itr_1 %>% filter(train_layer == 1) %>%  filter(test_layer == 2)
nonbinary_fidelity_merged_removed_itr_1_tn_1_ts_2 <- nonbinary_fidelity_merged_removed_itr_1_tn_1_ts_2 %>% filter(removed == 1)
nonbinary_fidelity_merged_removed_itr_1 <- nonbinary_fidelity_merged_removed_itr_1 %>% filter(removed == 1)
ggplot(nonbinary_fidelity_merged_removed_itr_1, aes(x = node_to, y = pred_obs)) +
  geom_boxplot() +
  theme_minimal() +
  coord_flip() +  # optional: flip if you have many species
  tme +
  labs(
    title = "Boxplot of pred_obs by pollinator species",
    x     = "Pollinator species (node_to)",
    y     = "pred_obs"
  )


## ---- check if the results relate to sample size ----

## ---- probability of observed zero interactions to actually be zero ----
df_all_itr <- read.csv('working_all_itr_nonbinary_canary.csv') # weighted
df_all_itr <- read.csv('working_all_itr_nonbinary_canary.csv') # predictions made on binary matrices

df_all_itr_zero <- df_all_itr %>% filter(removed == 1 & original_links == 0) # filter only removed links that are observed as zeros

# filter only links that appear several times for the analysis
df_summary <- df_all_itr_zero %>%
  group_by(node_from, node_to) %>%
  # Keep only those links/groups with >= 5 observations 
  filter(n() >= 5) %>%  
  summarise(
    n         = n(),
    mean_pred = mean(sigm_predicted, na.rm = TRUE),
    sd_pred   = sd(sigm_predicted, na.rm = TRUE),  # to check variance
    p_value   = if (sd_pred == 0) {
      NA_real_  # can't run a t-test if there's no variance
    } else {
      t.test(sigm_predicted, mu = 0.5)$p.value
    },
    .groups   = "drop"
  )

df_summary %>%
  arrange(desc(mean_pred)) %>%
  head(10)

write.csv(df_summary, 'links_probability_to_be_non_zeros.csv')

df_summary_sign <- df_summary %>% filter(p_value < 0.05)

df_top_20 <- df_summary %>%
  # Order by ascending p-value
  arrange(p_value) %>%
  # Take the first 20 rows
  slice(1:20)

ggplot(df_top_20, 
       aes(x = reorder(paste(node_to, node_from, sep = " - "), mean_pred),
           y = mean_pred,
           fill = p_value < 0.05)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_pred - sd_pred, ymax = mean_pred + sd_pred),
                width = 0.2) +
  coord_flip() +
  #scale_fill_manual(name = "Significant?", values = c("gray70", "tomato")) +
  theme_minimal() +
  scale_x_discrete(
    labels = function(x) sapply(x, function(lbl) {
      # 1) Replace underscores with a tilde (for spacing in plotmath)
      #    e.g., "Euphorbia_balsamifera_f" => "Euphorbia~balsamifera~f"
      lbl_tilde <- gsub("_", "~", lbl)
      # 2) Wrap in italic(), so the entire thing is in italics
      #    expression syntax: parse(text="italic(Euphorbia~balsamifera~f)")
      parse(text = paste0("italic(", lbl_tilde, ")"))
    })) +
  labs(
    x = "Link (pollinator - plant)",
    y = "Mean predicted value"
    #title = "Mean predicted value and significance"
  ) + tme

ggplot(df_top_20, aes(x = node_from, y = mean_pred)) +
  geom_col() +
  scale_x_discrete(
    labels = function(x) sapply(x, function(lbl) {
      # 1) Replace underscores with a tilde (for spacing in plotmath)
      #    e.g., "Euphorbia_balsamifera_f" => "Euphorbia~balsamifera~f"
      lbl_tilde <- gsub("_", "~", lbl)
      # 2) Wrap in italic(), so the entire thing is in italics
      #    expression syntax: parse(text="italic(Euphorbia~balsamifera~f)")
      parse(text = paste0("italic(", lbl_tilde, ")"))
    })
  ) +
  coord_flip() +
  theme_minimal()

## ---- train a model ----
set.seed(123)  # For reproducibility
n <- nrow(sites)
train_index <- sample(seq_len(n), size = floor(0.8 * n))  # 80% training
df_train <- sites[train_index, ]
df_test  <- sites[-train_index, ]

library(randomForest)

# Suppose your columns are named SVD_1, SVD_2, ..., SVD_k, and Weight is your target
# We'll use a simple formula interface for illustration:

rf_model <- randomForest(
  Weight ~ SVD_1 + SVD_2 + ... + SVD_k, 
  data = df_train,
  ntree = 500,      # Number of trees
  mtry = 3,         # Number of features tried at each split (tune this!)
  importance = TRUE # To get variable importance
)

