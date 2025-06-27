devtools::load_all()

# Define parameters
n <- 50                 # Total number of nodes
K <- 2                  # Number of communities
Z <- matrix(0, n, K)    # Membership matrix
Z[1:(n/2), 1] <- 1
Z[((n/2)+1):n, 2] <- 1

B <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)  # base connection probabilities
theta <- runif(n, 0.8, 1.2)  # degree heterogeneity

T <- 50                # Total time steps
t_star <- 30           # Time of structural change
persistence <- 0       # No change in community labels

# Time-varying B for local change in community 1
B_dynamic <- lapply(1:T, function(t) {
  B_curr <- B
  if (t >= t_star) {
    B_curr[1,1] <- 0.3  # local increase in interaction in community 1
  }
  return(B_curr)
})

# Generate dynamic SBM
dynamic_networks <- lapply(1:T, function(t) {
  generate_sbm(n=n, K=K, Z=Z, B=B_dynamic[[t]], theta=theta)
})


# Plot generation
# Pick time points
time_points <- c(10, 30, 45)

# Layout fixed across plots for consistent comparison
g_example <- graph_from_adjacency_matrix(dynamic_networks[[1]], mode = "undirected")
layout_pos <- layout_with_fr(g_example)

# Plot
par(mfrow = c(1, 3))  # 3 plots side by side
for (t in time_points) {
  g_t <- graph_from_adjacency_matrix(dynamic_networks[[t]], mode = "undirected")

  plot(g_t,
       layout = layout_pos,
       vertex.label = NA,
       vertex.size = 15,
       edge.color = "gray50",
       main = paste("t =", t))
}
par(mfrow = c(1, 1))  # reset
