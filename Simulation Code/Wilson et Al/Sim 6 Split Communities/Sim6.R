devtools::load_all()

# Define parameters
set.seed(42)
n <- 50                 # Total number of nodes
K <- 2                  # Number of communities
Z <- matrix(0, n, K)    # Membership matrix
Z[1:(n/2), 1] <- 1
Z[((n/2)+1):n, 2] <- 1

B <- matrix(c(0.2, 0.1, 0.1, 0.2), nrow=2, byrow=TRUE)  # base connection probabilities

theta <- runif(n, 0.8, 1.2)  # degree heterogeneity

T <- 50                # Total time steps
t_star <- 30           # Start of structural change
persistence <- 0       # No label change in this simulation


sim_output <- generate_dynamic_sbm(
  n = n,
  K = K,
  Z = Z,
  B = B,
  new_B = NULL,
  theta = theta,
  T = T,
  persistence = persistence,
  start_time = t_star,
  end_time = T,
  theta_fluctuate = FALSE,
  split_community = TRUE
)

dynamic_networks <- sim_output$adj_list
Z_list <- sim_output$Z_list

# Define label extractor
get_labels_from_Z <- function(Z) {
  apply(Z, 1, function(row) which(row == 1))
}

# Choose time points of interest
time_points <- c(10, 30, 45)

# Description for each time point
descriptions <- c(
  "t = 10:\nBefore structural change.",
  "t = 30:\nChange point activated.",
  "t = 45:\nPost-change."
)

# Fix layout for consistent visualization
#g_example <- graph_from_adjacency_matrix(dynamic_networks[[1]], mode = "undirected")
#layout_pos <- layout_with_fr(g_example)

# Manual layout: side-by-side positioning by community
initial_labels <- get_labels_from_Z(Z)
layout_pos <- matrix(NA, nrow = n, ncol = 2)
layout_pos[initial_labels == 1, 1] <- runif(sum(initial_labels == 1), min = -1, max = -0.2)
layout_pos[initial_labels == 2, 1] <- runif(sum(initial_labels == 2), min = 0.2, max = 1)
layout_pos[, 2] <- runif(n, min = -1, max = 1)  # vertical jitter


# Shewhart charts
par(mfrow = c(1, 3), mar = c(1, 1, 4, 1))  # More top margin for caption
for (i in 1:3) {
  t <- time_points[i]
  g_t <- graph_from_adjacency_matrix(dynamic_networks[[t]], mode = "undirected")
  Z_t <- Z_list[[t]]
  node_labels <- max.col((Z_t))

  color_map <- c("skyblue", "tomato", "gold", "gray")
  community_colors <- color_map[pmin(node_labels, length(color_map))]

  plot(g_t,
       layout = layout_pos,
       vertex.label = NA,
       vertex.size = 15,
       vertex.color = community_colors,
       edge.color = "gray50",
       main = descriptions[i])
}
par(mfrow = c(1, 1))

node_labels <- initial_labels

Z_list <- sim_output$Z_list

plot_simulation_summary(dynamic_networks, Z_list, node_labels, sim_title = "Simulation 6")

# === Compute F1–F9 summaries ===
F_time_series <- compute_all_F_series(dynamic_networks)

plot_F_summary_with_control_bands(F_time_series)
par(mfrow = c(1, 1))


lad_results <- run_lad_analysis(
  adjacency_list = dynamic_networks,
  changepoints = t_star,  # Known structural change
  title = "LAD - Simulation 6"
)
