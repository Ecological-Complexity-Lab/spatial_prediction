# ---- check k and lambda influence on f0.5 ----
# this pipeline allows us to predict missing links using the softImpute algorithm, calculate evaluators, have some stats and correlate the evaluators with ecological data.
# here we focus on island scale, but there is a section for comparison between scales.
# stages are according to the pipeline figure (Fig. 1).

# includes:
library(tidyverse)

source("code/common.R")


# 0) set an array of thresholds
thresholds <- seq(0, 10, by = 1)
thresholds <- thresholds/10 # create the doubles myself because seq function creates funky values for 0.6 and 0.7

# set lambda to filter
las <- c(1, 5, 50, 100)

# read data
alll <- readRDS("results/predictions_island_scale.rds")
all_ks <- alll %>% filter(input_lambda %in% las) %>%
  mutate(input_lambda = factor(input_lambda, levels = las))

# convert negatives to zeros
df <- all_ks %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# select the threshold for classifying links as 1s or 0s based on max F0.5
# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

# 2) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(train_layer, test_layer, k, itr, input_lambda, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f05_score   = (1.25) * (precision * recall) / ((0.25 * precision) + recall)
  ) %>%
  ungroup() 

# plot histograms of f05_score as a function of K and threshold
df_f05 <- df_thresh %>%
  select(k, itr,  input_lambda, threshold, f05_score)

# plot boxplot per K
df_f05 %>% filter(threshold == default_threshold) %>% 
ggplot(aes(x=as.factor(k), y=f05_score, color=as.factor(k))) +
  geom_boxplot() +
  facet_wrap(~input_lambda)





