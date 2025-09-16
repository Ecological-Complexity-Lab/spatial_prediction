# ---- partner fidelity affect on prediction ----
# this pipeline allows us to predict missing links of a single species across layers, 
# and so compare it (F1) to that species's partner fidelity

# Libraries -------------
library(tidyverse)
library(softImpute)
source("prediction_pipeline_for_publication/common.R") # load common functions

# Parameters -------------
k <- 2
prop_ones_to_remove <- 0.3 # proportion higher because the initial number is now (10 and up)
n_sim <- 50 # number of prediction iterations
degree_cutoff <- 7 # minimum degree to include a species in the analysis
n_layers_pf_cutoff <- 3 # minimum number of layers a species must be in to calculate partner fidelity
threshold <- 0.6 # threshold for converting probabilities to binary predictions
set.seed(42) # the answer to everything

# Themes ---------
tme <-  theme(axis.text = element_text(size = 18, color = "black"),
              axis.title = element_text(size = 18, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(color = "black"))
theme_set(theme_bw())

# Functions -------------
single_species_prediction <- function(networks, plant) {
  # function to predict links for a single species from one layer to another
  # returns a data frame with the results of all bootstrapping iterations
  
  # get list of layers this species has a degree of more then 10 in.
  populated_layers <- networks %>%
    filter(node_from == plant) %>%
    group_by(layer_from) %>%
    summarise(degree = n()) %>%
    filter(degree >= degree_cutoff) %>%
    pull(layer_from)
  
  combined_results <- data.frame()
  
  for (layers_to_train in populated_layers) {
    for (layer_to_predict in populated_layers) {
      if (layers_to_train == layer_to_predict) next # skip same layer prediction
      
      print(paste("** from:", layers_to_train, " to:", layer_to_predict, "**"))
      
      # Build the aggregated matrix A for training
      A <- build_interaction_matrix(data = networks, layers_to_filter = layers_to_train)
      
      # Build the layer to predict matrix P
      P <- build_interaction_matrix(data = networks, layers_to_filter = layer_to_predict)
      
      node_to <- rownames(P) # for the results
      node_from <- colnames(P)
                           
      ### ---- a. withhold links in P ----
      # map out the 0s and 1s in P - for this species only
      P_species <- P[, plant, drop = FALSE] # keep as matrix
      num_1_to_remove <- floor(sum(P_species>0, na.rm = T)*prop_ones_to_remove)  # Number of links to remove
      ones_in_P <- which(P_species > 0, arr.ind = TRUE)
      
      num_0_to_remove <- num_1_to_remove
      prop_0_removed <- num_0_to_remove / sum(P_species == 0, na.rm = T)
      zeros_in_P <- which(P_species == 0, arr.ind = TRUE)
      
      # make sure the inds are relevant to the whole matrix, not the sub-matrix
      plant_ind <- which(colnames(P) == plant)
      ones_in_P[,2] <- plant_ind
      zeros_in_P[,2] <- plant_ind
      
      # debug print
      print(paste("1 remove:", num_1_to_remove))
      print(paste("all 1   :", nrow(ones_in_P)))
      print(paste("0s to remove:", num_0_to_remove))
      print(paste("all zeros   :", nrow(zeros_in_P)))
      print(paste("prop of zeros removed   : ", prop_0_removed))
      
      # Randomly select zeros to withhold - bootstrapping
      bootstrapping_results <- NULL
      P_original <- P # save it for later
      
      for (i in 1:n_sim) {
        # remove 1s
        remove_indices <- ones_in_P[sample(1:nrow(ones_in_P), num_1_to_remove), ]
        P[remove_indices] <- NA  # Set removed links to NA
        
        # sample 0s
        zeros_to_remove_indices <- zeros_in_P[sample(1:nrow(zeros_in_P), num_0_to_remove), ]
        P[zeros_to_remove_indices] <- NA
        
        ### ---- creating a combined matrix C ----
        # Combine A and P into a single matrix C with NAs representing missing data
        all_row_ids <- unique(c(rownames(A), rownames(P)))
        all_col_ids <- unique(c(colnames(A), colnames(P)))
        C <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
                    dimnames = list(all_row_ids, all_col_ids))
        
        # Place A into C
        C[rownames(A), colnames(A)] <- A
        
        # Place P into C
        # Ensure that existing entries are not overwritten; sum overlapping entries
        C[rownames(P), colnames(P)] <- ifelse(is.na(C[rownames(P), colnames(P)]), 
                                              NA, 
                                              C[rownames(P), colnames(P)] + P[rownames(P), colnames(P)])
        
        
        # Apply biScale to center matrices
        C <- biScale(C, row.center=TRUE, col.center=TRUE, row.scale=FALSE, col.scale=FALSE)
        
        sum(is.na(C))
        
        ## ---- b. prediction with SVD ----
        k_values <- c(k)
        lam0 <- lambda0(C)
        lambda_values <- c(lam0)
        
        # Initialize variables to store the best results
        results <- data.frame(k = integer(),
                              lambda = numeric(),
                              original_links = numeric(),
                              predicted_values = numeric(),
                              input_lambda = numeric())
        not_removed_all <- NULL
        
        # Loop over all combinations of k and lambda
        for (k in k_values) {
          for (lambda in lambda_values) {
            # imputation
            r <- implement_impute(C, k, lambda, P, remove_indices, zeros_to_remove_indices, P_original)
            r$results$input_lambda <- lambda
            r$not_removed$input_lambda <- lambda
            results <- rbind(results, r$results)
            not_removed_all <- rbind(not_removed_all, r$not_removed)
          }
        }
        
        ### ---- save results for current k/lambda combination ----
        # After finishing the k/lambda loops, append the 'results' to 'combined_results'
        # ---- (D) Append to combined_results
        complete_edges_all <- rbind(results, not_removed_all)
        complete_edges_all$itr <- i
        bootstrapping_results <- rbind(bootstrapping_results, complete_edges_all)
        
        # reset P
        P <- P_original
      }
      
      combined_results <- rbind(
        combined_results,
        cbind(
          data.frame(
            emln_id = emln_id,
            train_layer = layers_to_train,
            test_layer = layer_to_predict,
            prop_ones_removed = prop_ones_to_remove,
            amount_of_removed_1 = num_1_to_remove,
            amount_of_removed_0 = num_0_to_remove,
            prop_0_removed = prop_0_removed
          ),
          bootstrapping_results
        )
      )
    }
  }
  return(combined_results)
}

calculate_PF_Sor <- function(x) {
  # check if number of layers is more then cutoff
  if (length(unique(x$layer_from)) < n_layers_pf_cutoff) {
    return(data.frame(layer1 = NA, layer2 = NA, sorensen = NA))
  }
  
  output <- NULL
  x_df <- as.data.frame(x)
  rownames(x_df) <- x$layer_from
  
  idx_pairs <- combn(x$layer_from, 2)
  for (pair_id in 1:ncol(idx_pairs)) {
    idx <- idx_pairs[, pair_id]
    l1 <- idx[1]
    l2 <- idx[2]
    
    p1 <- x_df[l1,"partners"][[1]]
    p2 <- x_df[l2,"partners"][[1]]
    a  <- length(intersect(p1, p2))          # shared partners
    b  <- length(setdiff(p1, p2))            # unique to first layer
    c  <- length(setdiff(p2, p1))            # unique to second layer
    sor_val <- 2 * a / (2 * a + b + c)
    
    output <- rbind(output, 
                    data.frame(layer1 = c(l1, l2), 
                               layer2 = c(l2, l1), 
                               sorensen = sor_val))
  }
  return(output)
}

# Load data -------------
# read networks
networks <- read.csv("prediction_pipeline_for_publication/results/network_island_scale.csv")

# check how many interlayer edges we have
networks %>% filter(layer_from != layer_to) %>% nrow() # 0 interlayer edges
# check how much data we will have after filtering
#View(networks %>%
#  group_by(node_to, layer_from) %>%
#  summarise(degree = n()))

# Run pipeline -------------

## calculate prediction - start with plants-------
plants <- unique(networks$node_from)

all_results <- NULL
for (plant in plants) {
  n_layers <- networks %>%
    filter(node_from == plant) %>%
    group_by(layer_from) %>%
    summarise(degree = n()) %>%
    filter(degree >= degree_cutoff) %>%
    pull(layer_from) %>% length()
  
  if (n_layers > 1) { # this analysis includes 7 plants
    print(paste0(plant, ": " , n_layers, " layers")) 
    ress <- single_species_prediction(networks, plant)
    
    ress$species <- plant
    ress <- ress %>% select(species, everything())
    all_results <- rbind(all_results, ress)
  }
  
}

# Add F1 score calculation

species_f1 <- all_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values)) %>% 
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  ) %>%
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(species, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
  ) %>%
  ungroup()

species_f1 %>% select(species, 
                      P_layer = test_layer, 
                      A_layer = train_layer, 
                      iter = itr, 
                      f1 = f1_score) %>% 
  write.csv("prediction_pipeline_for_publication/results/single_species_prediction_f1_scores.csv", row.names = FALSE)

f1_res <- species_f1 %>%
  group_by(species, train_layer, test_layer) %>%
  summarise(mean_f1 = mean(f1_score, na.rm = TRUE),
            sd_f1 = sd(f1_score, na.rm = TRUE)
  ) %>%
  ungroup()

f1_res

## check prediction on pollinators -------
polls <- unique(networks$node_to)

for (pol in polls) {
  n_layers <- networks %>%
    filter(node_to == pol) %>%
    group_by(layer_from) %>%
    summarise(degree = n()) %>%
    filter(degree >= degree_cutoff) %>%
    pull(layer_from) %>% length()
  
  if (n_layers > 1) { # this analysis includes 2 pollinators
    print(n_layers)
  }
} # only 2 pollinators survive the filter - not enough for analysis.


## ---- partner fidelity -------
df_fidelity <- networks %>% filter(weight != 0)

#### ---- calculate plants fidelity ----
# for each plant (node_from) and layer, gather the pollinators (node_to).
df_plant_partners <- df_fidelity %>%
  distinct(layer_from, node_from, node_to) %>%
  group_by(node_from, layer_from) %>%
  summarise(partners = list(unique(node_to)), .groups = "drop")

# For each plant, compute Sorensen similarity across all pairs of layers.
df_sorensen <- 
  df_plant_partners %>% 
  group_by(node_from) %>% 
  group_modify(~calculate_PF_Sor(.x))

# clean results
df_sorensen <- df_sorensen %>% drop_na() %>% 
  select(species = node_from , 
         train_layer = layer1, 
         test_layer = layer2, 
         sorensen)

# Plot results -------------
# combine prediction performance and partner fidelity
both <- f1_res %>%
  left_join(df_sorensen, by = c("species", "train_layer", "test_layer"))
 
# plot as scatter plot with trand line
ggplot(both, aes(x = sorensen, y = mean_f1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(
    title = "Relationship between Partner Fidelity and Link Prediction F1 Score",
    x = "Partner Fidelity (Sorensen Similarity)",
    y = "Mean F1 Score"
  ) +
  theme_minimal() + tme

# ---- pollinators with degree above 4 ----
degree_cutoff <- 4

all_results_polls <- NULL
for (pol in polls) {
  n_layers <- networks %>%
    filter(node_to == pol) %>%
    group_by(layer_from) %>%
    summarise(degree = n()) %>%
    filter(degree >= degree_cutoff) %>%
    pull(layer_from) %>% length()
  
  if (n_layers > 1) { # this analysis includes 9 pollinators
    print(paste0(pol, ": " , n_layers, " layers")) 
    ress_pol <- single_species_prediction(networks, pol)
    
    ress_pol$species <- pol
    ress_pol <- ress_pol %>% select(species, everything())
    all_results_polls <- rbind(all_results_polls, ress_pol)
  }
  
}

# Add F1 score calculation

species_f1_poll <- all_results_polls %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values)) %>% 
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  ) %>%
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by(species, train_layer, test_layer, itr) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
  ) %>%
  ungroup()

species_f1_poll %>% select(species, 
                      P_layer = test_layer, 
                      A_layer = train_layer, 
                      iter = itr, 
                      f1 = f1_score) %>% 
  write.csv("prediction_pipeline_for_publication/results/single_species_prediction_f1_scores_polls.csv", row.names = FALSE)

f1_res_poll <- species_f1_poll %>%
  group_by(species, train_layer, test_layer) %>%
  summarise(mean_f1 = mean(f1_score, na.rm = TRUE),
            sd_f1 = sd(f1_score, na.rm = TRUE)
  ) %>%
  ungroup()

f1_res_poll

# fidelity

df_fidelity <- networks %>% filter(weight != 0)

#### ---- calculate plants fidelity ----
# for each plant (node_from) and layer, gather the pollinators (node_to).
df_poll_partners <- df_fidelity %>%
  distinct(layer_from, node_from, node_to) %>%
  group_by(node_to, layer_from) %>%
  summarise(partners = list(unique(node_from)), .groups = "drop")

# For each plant, compute Sorensen similarity across all pairs of layers.
df_sorensen_polls <- 
  df_poll_partners %>% 
  group_by(node_to) %>% 
  group_modify(~calculate_PF_Sor(.x))

# clean results
df_sorensen_polls <- df_sorensen_polls %>% drop_na() %>% 
  select(species = node_to , 
         train_layer = layer1, 
         test_layer = layer2, 
         sorensen)

# Plot results -------------
# combine prediction performance and partner fidelity
both_poll <- f1_res_poll %>%
  left_join(df_sorensen_polls, by = c("species", "train_layer", "test_layer"))

# plot as scatter plot with trand line
ggplot(both_poll, aes(x = sorensen, y = mean_f1)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "rosybrown2", fill = "lightblue") +
  labs(
    title = "Pollinators",
    x = "Partner fidelity (Sorensen Similarity)",
    y = "Mean F1 score"
  ) +
  theme_minimal() + tme
