# add a null model with distance shuffling
# run shuffling according to what is written in the ipad


# load the distance matrix
dist_mat_kmm <- as.matrix(read.csv("results/distance_matrix_km.csv", row.names = 1))
colnames(dist_mat_kmm) <- rownames(dist_mat_kmm)
f05_mat <-  as.matrix(read.csv("results/f05_distance_matrix.csv", row.names = 1))
colnames(f05_mat) <- rownames(f05_mat)

# populate lower triangle with list of random values
shuff_dist_mat <- function(dist_mat) {
  suff_vals <- sample(as.dist(dist_mat))
  
  n_layers <- nrow(dist_mat)
  mm <- matrix(NA, nrow=n_layers, ncol=n_layers,
               dimnames=list(1:n_layers, 1:n_layers))
  rand_val_ind <- 0
  for(i in 1:nrow(mm)) for(j in 1:i) {
    if (i != j) {
      rand_val_ind <- rand_val_ind + 1
      mm[i,j] <- mm[j,i] <- suff_vals[rand_val_ind]
    }
  }
  as.dist(mm)
}

# test if they follow a normal distribution
shapiro.test(as.dist(f05_mat)) # p-value = 0.07614
shapiro.test(as.dist(dist_mat_km)) # p-value = 0.2949

p_val_dist <- c()
r2_val_dist <- c()
plot_points <- NA
# shuffle the distance matrix 500-1000 times
for (s in 1:500) {
  shuff_dist <- shuff_dist_mat(dist_mat_km)
  
  # calculate correlation with shuffled distance matrix
  corr_res <- cor.test(as.dist(f05_mat), as.dist(shuff_dist))
  
  # store the F0.5 value in a vector
  p_val_dist <- c(p_val_dist, corr_res$p.value)
  r2_val_dist <- c(r2_val_dist, corr_res$estimate^2)
  plot_points <- rbind(plot_points, data.frame(f05_val = as.dist(f05_mat), 
                                               dist_val = as.dist(shuff_dist),
                                               shuf_num = s))
}

plot_points$origin <- "shuffled"
plot_points <- rbind(plot_points, data.frame(f05_val = as.dist(f05_mat), 
                                             dist_val = as.dist(dist_mat_km),
                                             shuf_num = 0,
                                             origin = "observed"))

# plot sactter plot of F0.5 vs distance for shuffled and observed data
ggplot(plot_points, aes(x = dist_val, y = f05_val, color = origin)) +
  geom_point(alpha = 0.8) +
  labs(title = "Scatter plot of F0.5 vs Distance",
       x = "Distance", y = "F0.5 Score") +
  theme_minimal() +
  scale_color_manual(values = c("shuffled" = "blue", "observed" = "red"))


# calculate observed p-value
observed_corr <- cor.test(as.dist(f05_mat), as.dist(dist_mat_km))
observed_p_val <- observed_corr$p.value
#observed_r2_val <- observed_corr$estimate^2

# plot the distribution of p values with ggplot
p_val_df <- data.frame(p_value = p_val_dist)
ggplot(p_val_df, aes(x = p_value)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
  geom_vline(xintercept = observed_p_val, color = "red", linetype = "dashed") +
  labs(title = "Distribution of p-values from shuffled distance matrices",
       x = "p-value", y = "Frequency") +
  theme_minimal()

# calculate significant p value:
non_random_dist <- 
  sum(p_val_dist <= observed_p_val) / length(p_val_dist) # p-val = 0.012


