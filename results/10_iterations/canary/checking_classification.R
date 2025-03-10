correlation <- cor.test(nonbinary_removed$cut_weighted_predicted, nonbinary_removed$original_links, use = "complete.obs", method = "pearson")
correlation
# Extract correlation coefficient and p-value
r_value <- round(correlation$estimate, 3)
p_value <- formatC(correlation$p.value, format = "e", digits = 2)  # or round as you prefer
label_text <- paste0("r = ", r_value, ", p = ", p_value)

ggplot(nonbinary_removed, aes(x = original_links, y = cut_weighted_predicted)) +
  geom_point(color = "steelblue", alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "thistle") +
  labs(x = "Weight of original links", y = "Predicted value") +
  theme_minimal() +
  tme + 
  annotate("text",
           x = 70, y = 0.9,   # Adjust depending on your data range
           label = label_text,
           size = 4,
           color = "black") +
  theme(
    # Add a black frame around facet labels with thickness
    #strip.background = element_rect(color = "black", fill = "white", size = 1.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_line(color = "black")
  )


d <- read_csv('nonbinary_equal_0_1_removal_60_1.csv') # canary islands

d <- d %>%
  filter(removed == 1) %>% 
  filter(k == 2) %>% 
  filter(lambda == 0.1) %>%
  mutate(cut_weighted_predicted=case_when(predicted_values<0 ~ 0,
                                          predicted_values>1 ~ 1,
                                          TRUE ~ predicted_values)) %>% 
  mutate(predicted_bin = if_else(cut_weighted_predicted > 0.5, 1, 0))

result_summary <- d %>%
  group_by(emln_id, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_links == 1 & cut_weighted_predicted == 1),
    FN = sum(original_links == 1 & cut_weighted_predicted == 0),
    TN = sum(original_links == 0 & cut_weighted_predicted == 0),
    FP = sum(original_links == 0 & cut_weighted_predicted == 1),
    specificity = TN / (TN + FP),
    precision = TP / (TP + FP),
    recall = TP / (TP + FN),
    f1_score = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy = (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  ) %>%
  ungroup() %>%
  group_by(emln_id, train_layer, test_layer) %>%
  summarise(
    TP = TP,
    FN = FN,
    TN = TN,
    FP = FP,
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE)
  ) %>%
  ungroup()

net <- emln::load_emln(60) # canary islands
net$layers
net_name <- net$layers %>% select(layer_id, name)
net_name
net_name <- net_name %>%
  mutate(name = gsub("_", " ", name))

# Display the updated tibble
print(net_name)

# Perform left joins to replace IDs with names
result_summary <- result_summary %>%
  mutate(train_layer = as.integer(as.character(train_layer)),  # Ensure train_layer is an integer
         test_layer = as.integer(as.character(test_layer))) %>% # Ensure test_layer is an integer
  left_join(net_name, by = c("train_layer" = "layer_id")) %>% 
  rename(train_layer_name = name) %>%                          # Use a different name for clarity
  left_join(net_name, by = c("test_layer" = "layer_id")) %>%  
  rename(test_layer_name = name)                               # Use a different name for clarity

result_summary$diagonal <- result_summary$train_layer == result_summary$test_layer


# canaries 

layer_to_layer_plot_all_itr_canary <- 
  ggplot(result_summary, aes(x = train_layer_name, y = test_layer_name, fill = mcc)) +
  # First draw the entire heatmap with white borders for all tiles
  geom_tile(color = "white", linewidth = 0.1) +  
  # Then draw the diagonal tiles on top with black borders
  geom_tile(data = result_summary[result_summary$train_layer == result_summary$test_layer, ],
            color = "black", linewidth = 1.2) +  # Black borders only for diagonal tiles
  scale_fill_gradient(low = "steelblue2", high = "salmon2", na.value = "gray") +  # Set NA values to gray
  labs(x = "Training layer", y = "Predicted layer", fill = "mcc") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Minimize margins
    panel.background = element_blank(), #This ensures no panel background layers are drawn, which might add extra space.
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  # Rotate x-axis labels by 45 degrees
  ) +
  coord_fixed()

print(layer_to_layer_plot_all_itr_canary)


ggplot(d, aes(x=predicted_values))

test <- read.csv('aggregated_equal_0_1_removal_60_1_filtered.csv')
# test <- 
#   d %>% filter(original_links == 1) %>% filter(k == 2)

ggplot(test, aes(x=predicted_values))+
  geom_histogram(fill='lightsteelblue') + geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Island scale",
       x = "Predicted values",
       y = "Count") + tme
# 
# test %>% 
#   mutate(sig_pred = sigmoid(predicted_values)) %>% 
#   ggplot(aes(x=sig_pred)) + geom_histogram(fill = 'lightsteelblue')

test <- test %>% mutate(sigm_pred = sigmoid(predicted_values))

ggplot(test, aes(x=sigm_pred))+
  geom_histogram(fill='lightsteelblue') + geom_vline(xintercept = 0.53, linetype = "dashed") +
  labs(title = "Island scale",
       x = "Predicted values (logistic)",
       y = "Count") + tme

tanh_transform <- function(x){
  (tanh(x)+1)/2
}

clip_transform <- function(x){
  case_when(x<0~0,
                  x>1~1,
                  TRUE~x)
}

test %>% 
  mutate(tanh_pred=tanh_transform(predicted_values)) %>% 
  ggplot(aes(x=tanh_pred)) +
  geom_histogram(fill='thistle') +
  geom_vline(xintercept = 0.5) +
  labs(title = "Island scale",
       x = "Predicted values (tanh transformed)",
       y = "Count") + tme

test %>% 
  mutate(clip = clip_transform(predicted_values)) %>% 
  ggplot(aes(x=clip)) + 
  geom_histogram(fill='darkseagreen3') +
  labs(title = "Island scale",
       x = "Predicted values (>1)",
       y = "Count") + tme
