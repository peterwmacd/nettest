devtools::load_all()

# Define parameters
n <- 50                 # Total number of nodes
K <- 2                  # Number of communities
Z <- matrix(0, n, K)    # Membership matrix
Z[1:(n/2), 1] <- 1
Z[((n/2)+1):n, 2] <- 1

B <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)  # base connection probabilities
new_B <- B
new_B[1, 1] <- 0.3  # local increase in community 1

theta <- runif(n, 0.8, 1.2)  # degree heterogeneity

T <- 50                # Total time steps
t_star <- 30           # Start of structural change
persistence <- 0       # No label change in this simulation

# Generate dynamic SBM using your function
dynamic_networks <- generate_dynamic_sbm(
  n = n,
  K = K,
  Z = Z,
  B = B,
  new_B = new_B,
  theta = theta,
  T = T,
  persistence = persistence,
  start_time = t_star,
  end_time = T  # Maintain new_B through end
)
# Define label extractor
get_labels_from_Z <- function(Z) {
  apply(Z, 1, function(row) which(row == 1))
}
node_labels <- get_labels_from_Z(Z)

# Assign community colors
community_colors <- c("skyblue", "tomato")[node_labels]

# Choose time points of interest
time_points <- c(10, 30, 45)

# Description for each time point
descriptions <- c(
  "t = 10:\nBefore structural change.\nCommunities are balanced.",
  "t = 30:\nChange point activated.\nHigher density in community 1 (red).",
  "t = 45:\nPost-change.\nCommunity 1 remains more internally connected."
)

# Fix layout for consistent visualization
g_example <- graph_from_adjacency_matrix(dynamic_networks[[1]], mode = "undirected")
layout_pos <- layout_with_fr(g_example)

# Plot with captions
par(mfrow = c(1, 3), mar = c(1, 1, 4, 1))  # More top margin for caption
for (i in 1:3) {
  t <- time_points[i]
  g_t <- graph_from_adjacency_matrix(dynamic_networks[[t]], mode = "undirected")

  plot(g_t,
       layout = layout_pos,
       vertex.label = NA,
       vertex.size = 15,
       vertex.color = community_colors,
       edge.color = "gray50",
       main = descriptions[i])
}
par(mfrow = c(1, 1))
