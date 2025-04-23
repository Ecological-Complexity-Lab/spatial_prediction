# say your original data is in a matrix or data.frame `X`
# rows = species, columns = ecological variables

## ---- run ----
# run the algorithm with 1 iteration
# ---- softImpute for predicting removed links in empirical networks - all layer combinations ----
# including removal of zeros and ones. no bootstrapping yet
# in this code we remove zeros and ones and get the predicted values (not labels but values based on the svd from the softImpute fit function) in a dataframe for further exploration. We get results for a range of k and lambdas.
## ---- to do ----
#   CHECK IF COORDINATES IMPLY WHAT SHOULD BE INSTEAD OF NA - AND THEN WHAT IS THE MEANING OF REMOVING LINKS?
# check that there are no duplicates
# bootstrapping
# add aggregation (same format, different file)
# fix combined_results to include all k and lambda values for non-removed links
# ---- results so far -----
# results for all links and layer to layer predictions in combined_results_0.2_rem_values_nonbinary_all_edges2.csv.
## ---- load libraries ----
library(ggplot2)
library(emln)
library(pheatmap)
library(gridExtra)
library(dplyr)
library(ggrepel)
library(factoextra)
library(tidyverse)

## ---- functions ----
build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')
  
  # Step 3: Create the matrix with specific row and column species
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows
  
  # Initialize an empty matrix
  interaction_matrix <- matrix(0, nrow = length(species_to), ncol = length(species_from),
                               dimnames = list(species_to, species_from))
  
  # Populate the matrix with aggregated weights
  for (i in 1:nrow(aggregated_data)) {
    row <- aggregated_data$node_to[i]    # Rows represent 'node_to' species
    col <- aggregated_data$node_from[i]  # Columns represent 'node_from' species
    interaction_matrix[row, col] <- aggregated_data$weight[i]
  }
  
  return(interaction_matrix)
}

## ---- parameters ----
emln_id <- 60
layers_to_train <- 1
layer_to_predict <- 7

## ---- themes ----
tme <-  theme(axis.text = element_text(size = 14, color = "black"),
              axis.title = element_text(size = 14, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

## ---- run ----
### ---- load matrices ----

# Load matrices
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = TRUE, directed = FALSE)$layers_igraph
A_l <- d$extended

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(ggrepel)

# --- your existing preprocessing to get A_l and num_layers ---
A_l <- d$extended %>%
  mutate(
    layer_num = as.numeric(gsub("layer_", "", layer_from)),
    aggregated_layer = ifelse(
      layer_num %% 2 == 1,
      paste0("layer_", layer_num, "_", layer_num + 1),
      paste0("layer_", layer_num - 1, "_", layer_num)
    )
  )

#A_l <- A_l %>% filter(layer_from == layer_to)

# get the unique layers you actually want to loop over
num_layers <- length(unique(A_l$layer_from))

# get the vector of unique plant species
plants <- unique(A_l$node_from)

library(dplyr)

plants_df <- A_l %>% 
  distinct(node_from) %>% 
  arrange(node_from) %>% 
  rename(plant_species = node_from)

# to view as a tibble
print(plants_df)
write.csv(plants_df, 'plant_taxonomy.csv')

# if you just want it back as a character vector:
plants <- plants_df %>% pull(plant_species)
# initialize empty tibble to collect everything
combined_results <- NULL

# loop

for (layer_to_predict in 1:num_layers) {
  # build your P for this layer
  P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)
  
  # do the SVD
  svd_res     <- svd(P)
  poll_scores <- svd_res$u %*% diag(svd_res$d)
  plant_scores<- svd_res$v %*% diag(svd_res$d)
  
  # make two little tibbles
  poll_df <- tibble(
    species      = rownames(P),
    PC1          = poll_scores[,1],
    PC2          = poll_scores[,2],
    layer        = layer_to_predict,
    species_type = "pollinator"
  )
  
  plant_df <- tibble(
    species      = colnames(P),
    PC1          = plant_scores[,1],
    PC2          = plant_scores[,2],
    layer        = layer_to_predict,
    species_type = "plant"
  )
  
  # bind into the master table
  combined_results <- bind_rows(combined_results, poll_df, plant_df)
}

# At this point you have:
#   combined_results
#   # A tibble: N × 5
#      species   PC1    PC2    layer      species_type
#      <chr>    <dbl>  <dbl>  <chr>      <chr>       
#   1 spA      1.23    0.45   layer_1    pollinator  
#   2 spB     -0.67    2.11   layer_1    pollinator  
#   3 …        …       …      …          …           
#   …


# --- now manually create (or read in) a small lookup of families ---
# e.g.
# taxonomy_df <- tibble(
#   species = c("spA","spB","…"),
#   family  = c("Fam1","Fam2","…")
# )

## ---- add taxonomy ----
taxonomy_df <- read.csv('plant_taxonomy.csv')

combined_results <- combined_results %>%
  left_join(taxonomy_df, by = "species")
# now you have a ‘family’ column you can color by.


# --- example plotting ---
# 1) all species, colored by family, facetted by layer × species_type
ggplot(combined_results,
       aes(x = PC1, y = PC2, color = family, label = species)) +
  geom_point() +
  geom_text_repel(
    size = 2.5,
    max.overlaps = Inf,
    box.padding  = 0.3,
    point.padding= 0.2
  ) +
  facet_grid(species_type ~ layer) +
  theme_minimal() +
  labs(
    title = "PCA (SVD) ordinations by layer and species type",
    x     = "PC1",
    y     = "PC2"
  )

layer_1_plants <- combined_results %>% filter(species_type == "plant" & layer == "layer_1")
# 2) or, if you’d rather two separate series of small multiples:
## pollinators only
ggplot(filter(combined_results, species_type=="pollinator"),
       aes(PC1, PC2, color = family, label = species)) +
  geom_point() + geom_text_repel() +
  facet_wrap(~ layer) +
  theme_minimal() +
  labs(title = "Pollinator PCA by layer", x="PC1", y="PC2")


layer_1 <- combined_results %>% filter(layer == 1) # if we want to filter a certain layer
## plants only
ggplot(filter(layer_1, species_type=="plant"),
       aes(PC1, PC2, color = family, label = species)) +
  geom_point() + geom_text_repel() +
  facet_wrap(~ layer) +
  theme_minimal() +
  labs(title = "Plant PCA by layer", x="PC1", y="PC2") + tme
