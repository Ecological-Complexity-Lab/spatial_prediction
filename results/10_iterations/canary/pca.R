# say your original data is in a matrix or data.frame `X`
# rows = species, columns = ecological variables

res.pca <- prcomp(C_reconstructed,
                  center = TRUE,    # subtract variable means
                  scale.  = TRUE)   # divide by variable SDs

# Visualize with factoextra:
library(factoextra)
fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             pointsize  = 3,
             label      = "ind",    # or "none" if you only want points
             repel      = TRUE,
             title      = "PCA of pollinator species")

library(factoextra)
library(ggrepel)    # ggrepel is automatically used when you set repel = TRUE

# assume you've already done
# res.pca <- prcomp(X, center = TRUE, scale. = TRUE)

fviz_pca_ind(
  res.pca,
  geom.ind  = c("point", "text"),  # draw points *and* labels
  pointshape = 21,                 # filled circle
  pointsize  = 2,                  # size of the points
  label      = "all",              # label all “individuals” (i.e. your species)
  repel      = TRUE,               # use ggrepel to avoid text overlap
  labelsize  = 3,                  # font size for species names
  title      = "PCA of pollinator species"
)

# C_reconstructed: matrix (n_poll × n_plant)
svd_res <- svd(C_reconstructed)

# svd_res$u is n_poll × n_poll
# svd_res$v is n_plant × n_plant
# svd_res$d is length = min(n_poll,n_plant)

# Compute “scores” = U %*% Σ  and  V %*% Σ
poll_scores  <- svd_res$u %*% diag(svd_res$d)    # rows = pollinators
plant_scores <- svd_res$v %*% diag(svd_res$d)    # rows = plants

library(tibble)

poll_df <- as_tibble(
  poll_scores[,1:2],
  .name_repair = ~ c("PC1","PC2")
) %>%
  add_column(species = rownames(C_reconstructed), .before=1)

plant_df <- as_tibble(
  plant_scores[,1:2],
  .name_repair = ~ c("PC1","PC2")
) %>%
  add_column(species = colnames(C_reconstructed), .before=1)

library(ggplot2)
library(ggrepel)

# Pollinators
ggplot(poll_df, aes(x = PC1, y = PC2, label = species)) +
  geom_point() +
  geom_text_repel(size = 3) +
  labs(title = "SVD‐ordination: pollinators",
       x     = "Axis 1",
       y     = "Axis 2") +
  theme_minimal()

# Plants
ggplot(plant_df, aes(x = PC1, y = PC2, label = species)) +
  geom_point() +
  geom_text_repel(size = 3) +
  labs(title = "SVD‐ordination: plants",
       x     = "Axis 1",
       y     = "Axis 2") +
  theme_minimal()
