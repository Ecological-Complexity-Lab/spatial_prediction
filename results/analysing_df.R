# ---- exploring the df ----
# no bootstraping or aggregation yet
## ---- results so far ----
# binarization based on sigmoid is better than thresholding (0.5 == 0 otherwise 1).

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

### ---- trying to find the best way to convert the predicted values to weights ----
d_filtered <- d %<>% 
  mutate(binary_sigm=sigmoid(predicted_values)) %>% 
  mutate(binary_sigm=ifelse(binary_sigm>0.5,1,0)) %>% 
  mutate(binary_thresh=ifelse(predicted_values>0.5,1,0)) %>% 
  mutate(original_links_bin=ifelse(original_links>0.5,1,0)) %>% 
  mutate(TP=ifelse(binary_thresh==original_links_bin,1,0))
table(d$binary_sigm)
table(d$binary_thresh)

min_value <- min(d$predicted_values, na.rm = TRUE)

hist(d$predicted_values, 
     main = "Histogram of Predicted Values", 
     xlab = "Predicted Values", 
     ylab = "Frequency", 
     col = "blue", 
     border = "black")

table(d$TP)

d_removed <- filter(d, removed == 1)
d_removed <- d_removed %>%
  mutate(predicted_sigm = sigmoid(predicted_values)) %>%
  mutate(predicted_sigm_binary = ifelse(predicted_sigm>0.5, 1, 0)) %>% 
  mutate(original_links_binary = ifelse(original_links>=1, 1, 0)) %>% 
  mutate(TP = ifelse(predicted_sigm_binary == original_links_binary, 1, 0)) %>% 
  mutate(predicted_threshold = ifelse(predicted_sigm == 0.5, 0, 1)) %>% 
  mutate(theresh_TP = ifelse(predicted_threshold == original_links_binary, 1, 0))

table(d_removed$thresh_TP)

View(d_removed)

hist(d_removed$binary_sigm, 
     main = "Histogram of Predicted Values", 
     xlab = "Predicted Values", 
     ylab = "Frequency", 
     col = "blue", 
     border = "black")

# Pearson Correlation
pearson_cor <- cor(d_removed$predicted_values, d_removed$original_links, method = "pearson")

# Spearman Correlation
spearman_cor <- cor(d_removed$predicted_values, d_removed$original_links, method = "spearman")

# Create a data frame for plotting
df_plot <- data.frame(
  Observed = d_removed$original_links,  
  Predicted = d_removed$predicted_values
)

fit <- lm(Observed ~ Predicted, data = df_plot)
summary(fit)

alpha <- coef(fit)[1]          # Intercept
beta  <- coef(fit)[2]          # Coefficient for predicted_values

# New SVD-based predictions after fitting
df_plot$predicted_values_rescaled <- alpha + beta * df_plot$Predicted

# RMSE
rmse <- sqrt(mean((df_plot$Observed - df_plot$predicted_values_rescaled)^2))
cat("RMSE =", rmse, "\n")

# Correlation
cor(df_plot$Observed, df_plot$predicted_values_rescaled)



# Scatter plot with a line representing perfect prediction
ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Predicted vs. Observed Weights",
       x = "Observed Weights",
       y = "Predicted Weights") +
  theme_minimal()

mse <- mean((d_removed$predicted_values - d_removed$original_links)^2)

# Root Mean Squared Error (RMSE)
rmse <- sqrt(mse)

# Mean Absolute Error (MAE)
mae <- mean(abs(predicted_weights - observed_weights))

cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")

# ggplot(d, aes(predicted_values))+geom_histogram()

