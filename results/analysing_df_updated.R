# ---- exploring the results df ----
# no aggregation yet

## ---- results so far ----
# binarization based on thresholding (0.5 == 0 otherwise 1) is better than sigmoid. 7.42 trues:falses for thresholding, 0.7992887 trues:falses for sigmoid (we get more falses than trues). probably did it wrong.

## ---- to do ----
# find k and lambda for which the results are most accurate
# find the best way to binarize the predicted values
# average the predicted values for the same interactions 
## ---- load libraries ----
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pROC)

## ---- functions ----
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

## ---- load df ----
# d <- read_csv('combined_results_0.2_rem_values_nonbinary_all_edges_hpc.csv') # this file is really heavy so it takes a long time

# different format for input to reduce its size
filename = "combined_results_0.2_rem_values_nonbinary_all_edges_hpc.rds"
#saveRDS(d, file = filename, compress = TRUE)
#saveRDS(d %>% select(-node_to, -node_from), file = "test.rds",
#        compress = TRUE)

setwd("~/Library/CloudStorage/GoogleDrive-kesemoon@gmail.com/My Drive/חשב מסלול מחדש/Shai/research/svd_embeddings_lp/softimpute/hpc/results")

# d <- readRDS(filename)
# head(d)

## ---- working on one iteration cause the data is huge ----
# 
# d %>% 
#   filter(itr==1) %>% 
#   filter(k==2) %>% 
#   filter(lambda==0.01) %>% 
#   write_csv('test_df.csv')
# 
# # d %>% 
# #   filter(itr==1) %>% 
# #   filter(k==2) %>% 
# #   filter(lambda==0.1) %>% 
# #   filter(removed==1) %>% 
# #   ggplot(aes(original_links,predicted_values))+geom_point()
# # 
# # count(removed)
# 
# test_df <- read_csv('test_df.csv') # this is the one to use for all initial analysis
# 
# test_df %>% 
#   filter(removed==1) %>% 
#   ggplot(aes(original_links,predicted_values))+geom_point()
# 
# 
# ### ---- binarize the links ----
# 
# test_df_removed <- test_df %>% 
#   filter(removed==1) %>% 
#   write_csv('test_df_removed.csv')
# 
# test_df_removed <- test_df_removed %>%
#   mutate(original_bin = if_else(original_links > 0, 1, 0)) # binarize the original links
# 
# # now find the best way to binarize the predicted values
# # # Compute ROC using the raw predicted values against the binary truth
# # roc_obj <- roc(test_df$original_bin, test_df$predicted_values)
# # 
# # # Plot ROC curve
# # plot(roc_obj, main = "ROC Curve for Predicted Values")
# # 
# # # Calculate the optimal threshold based on Youden's J statistic (for example)
# # optimal_coords <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))
# # optimal_threshold <- optimal_coords["threshold"]
# # print(optimal_coords)
# # 
# # threshold0.5 <- 0.5
# # 
# # test_df <- test_df %>%
# #   mutate(predicted_bin_opt = if_else(predicted_values > optimal_threshold, 1, 0)) # binarize based on optimal threshold
# # 
# # test_df <- test_df %>%
# #   mutate(predicted_bin0.5 = if_else(predicted_values > threshold0.5, 1, 0)) # binarize based on threshold 0.5
# 
# # test_df <- test_df %>% 
# #   mutate(predicted_bin_eq_0.5 = if_else(predicted_values == threshold0.5, 0, 1))
# # 
# # table(test_df$original_bin, test_df$predicted_bin, dnn = c("Actual", "Predicted"))
# # 
# # table(test_df$original_bin, test_df$predicted_bin_opt, dnn = c("Actual", "Predicted"))
# # 
# # table(test_df$original_bin, test_df$predicted_bin0.5, dnn = c("Actual", "Predicted"))
# # 
# # table(test_df$original_bin, test_df$predicted_bin_eq_0.5, dnn = c("Actual", "Predicted"))
# 
# test_df_removed <- read_csv('test_df_removed.csv')
# 
# test_df_removed %>% 
#   ggplot(aes(original_links, predicted_prob)) + geom_point() # plotting the relationship between original links and predicted probabilities derived from logistic regression
# 
# # Basic Pearson correlation
# correlation_value <- cor(test_df_removed$predicted_prob, test_df_removed$original_links, method = "pearson")
# print(correlation_value)
# 
# spearman_value <- cor(test_df_removed$predicted_prob, test_df_removed$original_links, method = "spearman")
# print(spearman_value)
# 
# test_df_removed <- test_df_removed %>%
#   mutate(predicted_bin = if_else(predicted_prob > 0.5, 1, 0))
# 
# test_df_removed <- test_df_removed %>%
#   mutate(predicted_bin_0.3 = if_else(predicted_prob > 0.3, 1, 0)) # trying different thresholds
# 
# test_df_removed <- test_df_removed %>%
#   mutate(predicted_bin_0.8 = if_else(predicted_prob > 0.8, 1, 0))
# 
# table(test_df_removed$original_bin, test_df_removed$predicted_bin, dnn = c("Actual", "Predicted")) # we get a lot of false positives with the sigmoid function. 
# 
# table(test_df_removed$original_bin, test_df_removed$predicted_bin_0.3, dnn = c("Actual", "Predicted")) # reducing the threshold for binarization results in more positives (also false positives)
# 
# table(test_df_removed$original_bin, test_df_removed$predicted_bin_0.8, dnn = c("Actual", "Predicted")) # increasing the threshold results in more negatives (including false ones)
# 
# # trying min-max normalization instead of the logistic function
# min_val <- min(test_df_removed$predicted_values, na.rm = TRUE)
# max_val <- max(test_df_removed$predicted_values, na.rm = TRUE)
# 
# test_df_removed <- test_df_removed %>%
#   mutate(predicted_prob_minmax = (predicted_values - min_val) / (max_val - min_val))
# 
# test_df_removed %>% 
#   ggplot(aes(original_links, predicted_prob_minmax)) + geom_point()
# 
# spearman_value <- cor(test_df_removed$predicted_prob_minmax, test_df_removed$original_links, method = "spearman")
# print(spearman_value) # this is much worse than sigmoid. so sigmoid it is!
# 
# test_df_removed <- test_df_removed %>%
#   mutate(predicted_bin_minmax = if_else(predicted_prob_minmax > 0.5, 1, 0))
# 
# table(test_df_removed$original_bin, test_df_removed$predicted_bin_minmax, dnn = c("Actual", "Predicted")) # very bad!
# 
# ## ---- turn the original links into probabilities ----
# test_df_removed <- test_df_removed %>%
#   group_by(train_layer, test_layer, k, lambda, itr) %>%
#   mutate(relative_links = original_links / sum(original_links)) %>%
#   ungroup()

# here is a summarized table.

# test_df_removed <- test_df_removed %>% 
#   mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
#   mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) %>% 
#   mutate(original_bin = if_else(original_links > 0, 1, 0)) %>% # binarize the original links
#   group_by(train_layer, test_layer, k, lambda, itr) %>%
#   mutate(relative_links = original_links / sum(original_links)) %>%
#   ungroup() %>% 
#   write_csv('test_df_removed_binary_and_sigm.csv') # work with this table from here.

test_df_removed <- read_csv('combined_results_0.2_rem_values_nonbinary_all_edges_hpc.csv')

test_df_removed %>% 
  ggplot(aes(relative_links, predicted_values)) + geom_point() # plotting the relationship between original links and predicted probabilities derived from logistic regression

test_df_removed %>% 
  ggplot(aes(relative_links, predicted_prob_sigm)) + geom_point()

pearson <- cor(test_df_removed$predicted_values, test_df_removed$relative_links, method = "pearson") # very bad!

pearson <- cor(test_df_removed$predicted_prob_sigm, test_df_removed$relative_links, method = "pearson")

test_df_removed %>% 
  ggplot(aes(original_links, predicted_values)) + geom_point() # plotting the relationship between original links and predicted probabilities derived from logistic regression

## ---- some evaluators ----

TP <- sum(test_df_removed$original_bin == 1 & test_df_removed$predicted_bin == 1)
FN <- sum(test_df_removed$original_bin == 1 & test_df_removed$predicted_bin == 0)
TN <- sum(test_df_removed$original_bin == 0 & test_df_removed$predicted_bin == 0)
FP <- sum(test_df_removed$original_bin == 0 & test_df_removed$predicted_bin == 1)

specificity = TN / (TN + FP)
precision = TP / (TP + FP)
recall = TP / (TP + FN) # split to recall of removed and non-removed links
f1_score = 2 * (precision * recall) / (precision + recall)
balanced_accuracy = (recall + specificity) / 2
mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

evaluation_df <- data.frame(
  emln_id = emln_id,  
  train_layer = train_layer,  
  test_layer = test_layer,  
  k = k,
  lambda = lambda,
  TP = TP,
  TN = TN,
  FP = FP,
  FN = FN,
  specificity = specificity,
  precision = precision,
  recall = recall,
  f1_score = f1_score,
  balanced_accuracy = balanced_accuracy,
  mcc = mcc
)
# test_df %>% 
#   count(predicted_bin_opt) # using the threshold from ROC analysis gives all zeros...
# 
# test_df %>% 
#   count(predicted_bin)

# test_df %>% 
#   count(predicted_bin_eq_0.5)

## ---- now find the best 

# test_df %>% write_csv('test_df_0.5.csv')

## ---- get average probabilities for bootstrap iteration ----

## ---- start exploring ----

d_filtered_itr_1 <- d[d$itr == 1 & d$removed == 1 & d$original_links == 0,]

head(d_filtered_itr_1)
min(d_filtered_itr_1$predicted_values)
max(d_filtered_itr_1$predicted_values) # there is a large range of predicted values for the zeros in each iteration (-105.803 : 120.1398 for iteration 1)

d_filtered_itr_2 <- d[d$itr == 2 & d$removed == 1 & d$original_links == 0,]

head(d_filtered_itr_2)

identical(d_filtered_itr_1$node_to, d_filtered_itr_2$node_to)

d_filtered_itr_1_ones <- d[d$itr == 1 & d$removed == 1 & d$original_links != 0,]
d_filtered_itr_2_ones <- d[d$itr == 2 & d$removed == 1 & d$original_links != 0,]
identical(d_filtered_itr_1_ones$node_to, d_filtered_itr_2_ones$node_to) # confirming that in each iteration only the zeros are resampled.

### ---- trying to find the best way to convert the predicted values to weights ----
d_filtered <- d %<>% 
  mutate(binary_sigm=sigmoid(predicted_values)) %>% 
  mutate(binary_sigm=ifelse(binary_sigm>0.5,1,0)) %>% 
  mutate(binary_thresh=ifelse(predicted_values>0.5,1,0)) %>% 
  mutate(original_links_bin=ifelse(original_links>0.5,1,0)) %>% 
  mutate(TP_thresh=ifelse(binary_thresh==original_links_bin,1,0)) %>%
  mutate(TP_sigm=ifelse(binary_sigm==original_links_bin,1,0))

table(d$TP_thresh)
table(d$TP_sigm)

d_1_1 <- d[d$train_layer == 1 & d$test_layer == 1,]
fit <- lm(original_links ~ predicted_values, data = d_1_1)
summary(fit)


table(d$binary_sigm)
table(d$binary_thresh)

min_value <- min(d$predicted_values, na.rm = TRUE)

hist(d$predicted_values, 
     main = "Histogram of Predicted Values", 
     xlab = "Predicted Values", 
     ylab = "Frequency", 
     col = "blue", 
     border = "black")

table(d$TP_thresh)
table(d$TP_sigm)

d_removed <- filter(d, removed == 1)
d_removed <- d_removed %>%
  mutate(predicted_sigm = sigmoid(predicted_values)) %>%
  mutate(predicted_sigm_binary = ifelse(predicted_sigm>0.5, 1, 0)) %>% 
  mutate(original_links_binary = ifelse(original_links>=1, 1, 0)) %>% 
  mutate(TP = ifelse(predicted_sigm_binary == original_links_binary, 1, 0)) %>% 
  mutate(predicted_threshold = ifelse(predicted_sigm == 0.5, 0, 1)) %>% 
  mutate(theresh_TP = ifelse(predicted_threshold == original_links_binary, 1, 0))

table(d_removed$thresh_TP)

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

## ---- Shai ----
# Define a function that computes evaluation metrics from a dataframe with binary obs and pred
eval_metrics <- function(data) {
  data %>%
    summarise(
      TP = sum(obs == 1 & pred == 1),
      FP = sum(obs == 0 & pred == 1),
      FN = sum(obs == 1 & pred == 0),
      TN = sum(obs == 0 & pred == 0)
    ) %>%
    mutate(
      precision         = ifelse((TP + FP) > 0, TP / (TP + FP), NA),
      recall            = ifelse((TP + FN) > 0, TP / (TP + FN), NA),
      f1                = ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                                 2 * precision * recall / (precision + recall), NA),
      true_negative_rate = ifelse((TN + FP) > 0, TN / (TN + FP), NA),
      balanced_accuracy = (recall + true_negative_rate) / 2
    )
}

test_df_removed <- 
  d %>%
  filter(k==5) %>%
  filter(lambda==0.01) %>%
  filter(removed==1)



test_df_removed %>% 
  # filter(original_links>0) %>% 
  mutate(obs=ifelse(original_links>0,1,0)) %>% 
  mutate(predicted_prob_sigm = sigmoid(predicted_values)) %>%  # convert the predicted values to probability values in the interval (0, 1) using the logistic function
  mutate(predicted_bin_sigm = if_else(predicted_prob_sigm > 0.5, 1, 0)) %>%
  
  # mutate(pred=ifelse(predicted_values>0,1,0)) %>% 
  mutate(pred=predicted_bin_sigm) %>% 
  group_by(itr,train_layer, test_layer) %>% 
  eval_metrics() %>% 
  ungroup() %>% 
  group_by(train_layer, test_layer) %>% 
  summarise(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    mean_f1 = mean(f1, na.rm = TRUE),
    mean_balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE)
  ) %>% 
  
  ggplot(aes(x=as.factor(train_layer), y=as.factor(test_layer), fill=mean_f1))+
  geom_tile()+
  scale_fill_gradient2(
    low = "red",      # Color for low values
    mid = "white",    # Color at the midpoint
    high = "blue",    # Color for high values
    midpoint = 0.5    # Set white at 0.5
  )


test_df_removed %>% 
  distinct(test_layer, train_layer)


