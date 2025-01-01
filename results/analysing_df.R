# ---- exploring the df ----
# no bootstrapin or aggregation yet

## ---- load libraries ----
library(tidyverse)
library(ggplot2)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load df ----
d <- read_csv('combined_results_0.2_rem_values_nonbinary_all_edges2.csv')
d
View(d)

## ---- start exploring ----

d %<>% 
  mutate(binary_sigm=sigmoid(predicted_values)) %>% 
  mutate(binary_sigm=ifelse(binary_sigm>0.5,1,0)) %>% 
  mutate(binary_thresh=ifelse(predicted_values>0.5,1,0)) %>% 
  mutate(original_links_bin=ifelse(original_links>0.5,1,0)) %>% 
  mutate(TP=ifelse(binary_thresh==original_links_bin,1,0))
table(d$binary_sigm)
table(d$binary_thresh)

length(d$predicted_values)

min_value <- min(d$predicted_values, na.rm = TRUE)

hist(d$predicted_values, 
     main = "Histogram of Predicted Values", 
     xlab = "Predicted Values", 
     ylab = "Frequency", 
     col = "blue", 
     border = "black")

table(d$TP)

# Pearson Correlation
pearson_cor_all <- cor(predicted_all, observed_all, method = "pearson")

# Spearman Correlation
spearman_cor_all <- cor(predicted_all, observed_all, method = "spearman")

cat("Pearson Correlation (All):", pearson_cor_all, "\n")
cat("Spearman Correlation (All):", spearman_cor_all, "\n")

# Create a data frame for plotting
df_plot <- data.frame(
  Observed = observed_weights,  
  Predicted = predicted_weights
)

mse <- mean((predicted_weights - observed_weights)^2)

# Root Mean Squared Error (RMSE)
rmse <- sqrt(mse)

# Mean Absolute Error (MAE)
mae <- mean(abs(predicted_weights - observed_weights))

cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")

# ggplot(d, aes(predicted_values))+geom_histogram()

