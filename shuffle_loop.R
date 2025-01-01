# ---- compare consistency of svd feature vectors of pollinators ----
## in this code we calculate the Euclidean distance between a feature vector of a species of pollinator in a network layer to its vector in another layer, and compare it to the distance to the vectors of the same species in shuffled networks
## ---- load libraries ----
# install.packages("softImpute")
library(softImpute)
library(PRROC)
library(pROC)
library(ggplot2)
library(emln)
library(pheatmap)

## ---- functions ----

build_interaction_matrix <- function(data, layers_to_train) {
  
  # Step 1: Filter rows based on specified layers in 'layer_from'
  layers_to_filter <- paste0("layer_", layers_to_train)
  filtered_data <- subset(data, layer_from %in% layers_to_filter)
  
  # Step 2: Aggregate weights for identical species pairs
  library(dplyr)
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = 'drop')
  
  # Step 3: Create the matrix with specific row and column species
  # Define unique species in node_from and node_to for matrix dimensions
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to <- unique(aggregated_data$node_to)      # Rows
  
  # Initialize an empty matrix with rows as node_to species and columns as node_from species
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

join_matrices_non_overlapping <- function(A, B) {
  # Determine the size of the resulting matrix
  rows_A <- nrow(A)
  cols_A <- ncol(A)
  rows_B <- nrow(B)
  cols_B <- ncol(B)
  
  # Create a zero matrix with appropriate size
  result_matrix <- Matrix(0, nrow = rows_A + rows_B, ncol = cols_A + cols_B)
  
  # Place matrix A in the top-left corner
  result_matrix[1:rows_A, 1:cols_A] <- A
  
  # Place matrix B in the bottom-right corner
  result_matrix[(rows_A + 1):(rows_A + rows_B), (cols_A + 1):(cols_A + cols_B)] <- B
  
  # Combine row names
  rownames(result_matrix) <- c(rownames(A), rownames(B))
  
  # Combine column names
  colnames(result_matrix) <- c(colnames(A), colnames(B))
  
  return(result_matrix)
}

## ---- parameters ----
emln_id <- 25

## ---- run ----
### ---- load matrices ----
d <- load_emln(emln_id)
graph_list <- get_igraph(d, bipartite = T, directed = F)$layers_igraph # This includes in each layer all the nodes in the system (so could be many singletons)
# Total number of layers
num_layers <- length(graph_list)
A_l <- d$extended

# Initialize a list to store species from each network
species_list <- list()

# Collect species names from each network
for (layers_to_train in 1:num_layers) {
  P_temp <- build_interaction_matrix(data = A_l, layers_to_train = layers_to_train)
  P_temp[P_temp > 0] <- 1  # Make binary
  species_list[[layers_to_train]] <- rownames(P_temp)  # Assuming pollinators are in rows
}

# Find species common to all networks
common_species <- Reduce(intersect, species_list)

# Initialize an empty data frame to store all Euclidean distances
distances_all <- data.frame()

# Loop through all combinations of layers_to_train and layers_to_compare
for (layers_to_train in 1:num_layers) {
  for (layers_to_compare in 1:num_layers) {
    # Skip if comparing the same network
    if (layers_to_compare == layers_to_train) next
    
    # Build interaction matrices for both networks
    P <- build_interaction_matrix(data = A_l, layers_to_train = layers_to_train)
    P[P > 0] <- 1  # Make binary
    
    P_compare <- build_interaction_matrix(data = A_l, layers_to_train = layers_to_compare)
    P_compare[P_compare > 0] <- 1  # Make binary
    
    # Filter matrices to include only overlapping species
    P <- P[common_species, , drop = FALSE]
    P_compare <- P_compare[common_species, , drop = FALSE]
    
    # Shuffle the interactions while keeping species names
    P_shuffled <- P
    P_shuffled[] <- sample(P)
    
    P_compare_shuffled <- P_compare
    P_compare_shuffled[] <- sample(P_compare)
    
    ## ---- SVD ----
    # Combine P and P_compare to ensure a common embedding space
    AP <- join_matrices_non_overlapping(A = P, B = t(P_compare))
    
    # Perform SVD with k = 10
    k <- 10
    svd_AP <- svd(AP)
    
    svd_P <- svd_AP
    svd_P$u <- svd_P$u[1:nrow(P),]
    svd_P$v <- svd_AP$v[1:ncol(P), ]
    
    svd_P_shuffled <- svd_AP
    svd_P_shuffled$u <- svd_P_shuffled$u[(nrow(P_shuffled)+1):nrow(svd_P$u),]
    svd_P_shuffled$v <- svd_AP$v[(ncol(P_shuffled) + 1):ncol(AP), ]
    
    # Step 4: Reduce the embeddings to k dimensions for more efficient modeling and avoiding overfitting.
    # We retain only the top k dimensions, as selected earlier, to reduce complexity while preserving structure.
    svd_P_shuffled$u <- svd_P_shuffled$u[, 1:k]
    svd_P_shuffled$v <- svd_P_shuffled$v[, 1:k]
    svd_P_shuffled$d <- svd_P_shuffled$d[1:k]
    
    # Set row and column names
    rownames(svd_AP$u) <- rownames(AP)
    rownames(svd_AP$v) <- colnames(AP)
    
    # Extract feature vectors for P and P_compare
    num_rows_P <- nrow(P)
    num_cols_P <- ncol(P)
    num_rows_P_compare <- nrow(P_compare)
    num_cols_P_compare <- ncol(P_compare)
    
    svd_P <- list(
      u = svd_AP$u[1:num_rows_P, ],  # Pollinator features for P
      v = svd_AP$v[1:num_cols_P, ]   # Plant features for P
    )
    
    svd_P_compare <- list(
      u = svd_AP$u[(num_rows_P + 1):(num_rows_P + num_rows_P_compare), ],
      v = svd_AP$v[(num_cols_P + 1):(num_cols_P + num_cols_P_compare), ]
    )
    
    # Repeat for shuffled matrices
    AP_shuffled <- join_matrices_non_overlapping(A = P_shuffled, B = t(P_compare_shuffled))
    svd_AP_shuffled <- sparsesvd(AP_shuffled, rank = k)
    
    rownames(svd_AP_shuffled$u) <- rownames(AP_shuffled)
    rownames(svd_AP_shuffled$v) <- colnames(AP_shuffled)
    
    svd_P_shuffled <- list(
      u = svd_AP_shuffled$u[1:num_rows_P, ],
      v = svd_AP_shuffled$v[1:num_cols_P, ]
    )
    
    svd_P_compare_shuffled <- list(
      u = svd_AP_shuffled$u[(num_rows_P + 1):(num_rows_P + num_rows_P_compare), ],
      v = svd_AP_shuffled$v[(num_cols_P + 1):(num_cols_P + num_cols_P_compare), ]
    )
    
    ## ---- Euclidean distances ----
    # Extract pollinator feature vectors
    features_P <- svd_P$u
    features_P_compare <- svd_P_compare$u
    
    rownames(features_P) <- rownames(P)
    rownames(features_P_compare) <- rownames(P_compare)
    
    # Species names should match since we filtered for common_species
    for (species in common_species) {
      vec_P <- features_P[species, ]
      vec_P_compare <- features_P_compare[species, ]
      
      # Calculate Euclidean distance
      euclidean_distance <- sqrt(sum((vec_P - vec_P_compare)^2))
      
      # Store the results
      distances_all <- rbind(distances_all, data.frame(
        Species = species,
        Network1 = layers_to_train,
        Network2 = layers_to_compare,
        Euclidean_Distance = euclidean_distance,
        Shuffled = FALSE
      ))
    }
    
    # Repeat for shuffled matrices
    features_P_shuffled <- svd_P_shuffled$u
    features_P_compare_shuffled <- svd_P_compare_shuffled$u
    
    rownames(features_P_shuffled) <- rownames(P_shuffled)
    rownames(features_P_compare_shuffled) <- rownames(P_compare_shuffled)
    
    for (species in common_species) {
      vec_P_shuffled <- features_P_shuffled[species, ]
      vec_P_compare_shuffled <- features_P_compare_shuffled[species, ]
      
      # Calculate Euclidean distance
      euclidean_distance_shuffled <- sqrt(sum((vec_P_shuffled - vec_P_compare_shuffled)^2))
      
      # Store the results
      distances_all <- rbind(distances_all, data.frame(
        Species = species,
        Network1 = layers_to_train,
        Network2 = layers_to_compare,
        Euclidean_Distance = euclidean_distance_shuffled,
        Shuffled = TRUE
      ))
    }
  }
}
# Statistical Analysis
# Separate data for P and P_shuffled
distances_P <- subset(distances_all, Shuffled == FALSE)
distances_P_shuffled <- subset(distances_all, Shuffled == TRUE)

# Calculate means and variances
mean_P <- mean(distances_P$Euclidean_Distance)
var_P <- var(distances_P$Euclidean_Distance)

mean_P_shuffled <- mean(distances_P_shuffled$Euclidean_Distance)
var_P_shuffled <- var(distances_P_shuffled$Euclidean_Distance)

cat("Mean Euclidean Distance for P:", mean_P, "\n")
cat("Variance of Euclidean Distance for P:", var_P, "\n\n")

cat("Mean Euclidean Distance for P_shuffled:", mean_P_shuffled, "\n")
cat("Variance of Euclidean Distance for P_shuffled:", var_P_shuffled, "\n")

# Statistical tests
# Check normality
shapiro_P <- shapiro.test(distances_P$Euclidean_Distance)
shapiro_P_shuffled <- shapiro.test(distances_P_shuffled$Euclidean_Distance)

# Choose appropriate test based on normality
if (shapiro_P$p.value > 0.05 && shapiro_P_shuffled$p.value > 0.05) {
  # Use t-test
  t_test_result <- t.test(distances_P$Euclidean_Distance, distances_P_shuffled$Euclidean_Distance)
  print(t_test_result)
} else {
  # Use Wilcoxon test
  wilcox_test_result <- wilcox.test(distances_P$Euclidean_Distance, distances_P_shuffled$Euclidean_Distance)
  print(wilcox_test_result)
}

# Test for variance difference
levene_test_result <- leveneTest(Euclidean_Distance ~ Shuffled, data = distances_all)
print(levene_test_result)
