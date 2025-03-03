# ---- variable contribution ----
# load libaries
library(corrplot)

# load data
df <- read.csv("working_df_island_distance_fidelity_jaccard_netsize.csv")

# analyze only one off-diagonal
df_off <- df %>%
  # Keep rows where train_layer < test_layer (upper triangle) or on the diagonal
  filter(train_layer < test_layer)
# subset relavant columns
df_subset <- df_off %>% select(balanced_accuracy, distance_km,	avg_sorensen_plants,	avg_sorensen_pollinators,	jaccard_pollinators,	jaccard_plants,	jaccard_edges, size_P,	density_P, size_C,	density_C)

## ---- autocorrelation check ----
# 2. Compute correlation matrix
cor_mat <- cor(df_subset, use = "complete.obs")
#write.csv(cor_mat, "island_autocorrelation.csv")

# 3. Visualize (optional)

corrplot(cor_mat, 
         method = "number",      # or "circle", "color", etc.
         type = "upper",         # upper/lower/full
         tl.cex = 0.7,           # text label size
         number.cex = 0.7,       # correlation coefficient size
         tl.col = "black"       # text label color (optional)
)

## ---- variable importance ----

### ---- linear regression ----
# Fit a linear model predicting f1_score from all other numeric predictors
lm_fit <- lm(f1_score ~ distance_km +	avg_sorensen_plants +	avg_sorensen_pollinators +	jaccard_pollinators +	jaccard_plants +	jaccard_edges + size_P +	density_P + size_C +	density_C,
             data = df_subset)

summary(lm_fit)

### ---- random forest ----
library(randomForest)

# Random Forest (only with numeric columns)
rf_fit <- randomForest(balanced_accuracy ~ ., data = df_subset, importance = TRUE)

# Check variable importance
importance(rf_fit)
varImpPlot(rf_fit)

# Extract importance
imp <- importance(rf_fit) 
# For regression: imp is a matrix with columns: %IncMSE, IncNodePurity
# For classification: imp often has two columns per measure.
# We'll assume %IncMSE and IncNodePurity are present.

# Turn it into a data frame for easier plotting
imp_df <- as.data.frame(imp)
imp_df$F1_score <- rownames(imp_df)  # Keep variable names in a column

# Example for %IncMSE
ggplot(imp_df, aes(x = reorder(F1_score, `%IncMSE`), y = `%IncMSE`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Balanced accuracy", 
       y = "% Increase in MSE",
       title = "Variable Importance (by Permutation)") +
  tme

### ---- pca ----
# Install if needed
library(factoextra)

df_subset_explan <- df_subset %>% select(-f1_score)

pca_res <- prcomp(df_subset, scale. = TRUE)

# Visualize individuals (rows in your dataset)
fviz_pca_ind(
  pca_res,
  # color points by their cos2 (quality of representation)
  col.ind = "cos2", 
  # or color by a grouping variable in your original data frame
  # habillage = df$someFactor,
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE    # Avoid overlapping text labels
)

# Visualize variables (columns / features)
fviz_pca_var(
  pca_res, 
  col.var = "contrib",  # Color by contributions to the PC
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
) + tme

# Biplot (both individuals + variables)
fviz_pca_biplot(
  pca_res,
  repel = TRUE,
  col.var = "contrib", 
  col.ind = "cos2"
) +tme



