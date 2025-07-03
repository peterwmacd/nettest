# --- Static Kidney-Egg Plot ---
cat("Generating static Kidney-Egg graph...\n")

n <- 20
m <- 5
p <- 0.1
q <- 0.5
theta <- runif(n, 0.8, 1.2)
seed <- 123

# Generate static graph and extract true egg nodes
set.seed(seed)
egg_nodes <- sample(1:n, m)
kidney_nodes <- setdiff(1:n, egg_nodes)

A_static <- generate_kidney_egg(n, m, p, q, theta = theta, seed = seed)
g_static <- graph_from_adjacency_matrix(A_static, mode = "undirected")

# Color and label
node_colors <- rep("skyblue", n)
node_colors[egg_nodes] <- "tomato"
node_labels <- ifelse(1:n %in% egg_nodes, "E", "K")

layout_static <- layout_with_fr(g_static)

plot(g_static,
     layout = layout_static,
     main = "Static Kidney-Egg Graph",
     vertex.color = node_colors,
     vertex.label = node_labels,
     vertex.size = 20,
     vertex.label.color = "black",
     edge.color = "gray50")

# --- Dynamic Kidney-Egg Sequence ---
cat("Generating dynamic Kidney-Egg graph sequence...\n")

T <- 10
anomaly_start <- 4
anomaly_end <- 6

A_list <- generate_dynamic_kidney_egg(n, m, p, q, T,
                                      anomaly_start, anomaly_end,
                                      theta = theta, seed = seed)

# Use consistent layout
layout_dyn <- layout_static  # reuse for visual consistency

# Plot before, during, after
par(mfrow = c(1, 3))

time_points <- c(3, 5, 8)
titles <- c(
  "Before Anomaly (all background)",
  "During Anomaly (egg chatter active)",
  "After Anomaly (all background)"
)

par(mfrow = c(1, 3))
for (i in seq_along(time_points)) {
  t <- time_points[i]
  g_t <- graph_from_adjacency_matrix(A_list[[t]], mode = "undirected")

  # Color and label logic based on time
  if (t >= anomaly_start && t <= anomaly_end) {
    vcolor <- rep("skyblue", n)
    vcolor[egg_nodes] <- "tomato"
    vlabel <- ifelse(1:n %in% egg_nodes, "E", "K")
  } else {
    vcolor <- rep("skyblue", n)
    vlabel <- rep("K", n)
  }

  plot(g_t,
       layout = layout_dyn,
       main = titles[i],
       vertex.color = vcolor,
       vertex.label = vlabel,
       vertex.size = 20,
       vertex.label.color = "black",
       edge.color = "gray60")
}
par(mfrow = c(1, 1))

