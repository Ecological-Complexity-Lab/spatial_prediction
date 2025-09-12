A_l_34 <- A_l %>% filter(layer_from == "layer_3" | layer_from == "layer_4")
agg_4 <- aggregated_df %>% filter(layer_from == "layer_4")
unique(A_l_34$node_from)
unique(agg_4$node_from)
intersect(A_l_34$node_from, agg_4$node_from)
unique(A_l_34$node_to)
unique(agg_4$node_to)
intersect(A_l_34$node_to, agg_4$node_to)

data <- A_l
layers <- "layer_1"
subset(data, layer_from %in% layers)

A_l12 <- A_l %>% filter(layer_from == "layer_1" | layer_from == "layer_2")
agg_1 <- aggregated_df %>% filter(layer_from == "layer_1")
unique(A_l12$node_from)
unique(agg_1$node_from)
intersect(A_l12$node_from,agg_1$node_from)

A_l34 <- A_l %>% filter(layer_from == "layer_3" | layer_from == "layer_4")
agg_2 <- aggregated_df %>% filter(layer_from == "layer_2")
unique(A_l34$node_from)
unique(agg_2$node_from)
intersect(A_l34$node_from,agg_2$node_from) # the species in aggregated df are aligned with the ones in corresponding layers in A_l

# now i ran the pipeline until the prediction loop
layer_to_predict <- 1
# Build the aggregated matrix A for training
A <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layers_to_train)
# Build the layer to predict matrix P
P <- build_interaction_matrix(data = aggregated_df, layers_to_filter = layer_to_predict)
node_to <- rownames(P) # for the results
node_from <- colnames(P)
unique(node_from)
intersect(agg_1$node_from,node_from) # all species are there
unique(node_to)
intersect(agg_1$node_to,node_to) # same.

# now check site scale
# i ran the pipeline until the prediction loop
layer_to_predict <- 1
layers_to_train <- 7
# Build the aggregated matrix A for training
A <- build_interaction_matrix(data = A_l, layers_to_filter = layers_to_train)

# Build the layer to predict matrix P
P <- build_interaction_matrix(data = A_l, layers_to_filter = layer_to_predict)

node_to <- rownames(P) # for the results
node_from <- colnames(P)
unique(node_from)
A_l1 <- A_l %>% filter(layer_from == "layer_1")
unique(A_l1$node_from)
intersect(A_l1$node_from,node_from) # species are the same.

# check if the dimensions of P match the ones in ST1
# P is layer 7, which is La Gomera according to d$layers
agg7 <- aggregated_df %>% filter(layer_from == "layer_7")
length(unique(agg7$node_from))
length(unique(agg7$node_to)) # we have 67 species which is the same as in ST1 yay
