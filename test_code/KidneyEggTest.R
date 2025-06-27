devtools::load_all()


# ----- Test: Static Kidney-Egg Graph -----
cat("Generating static Kidney-Egg graph...\n")

n <- 20         # total nodes
m <- 5          # egg (chatter) nodes
p <- 0.1        # background (kidney) edge probability
q <- 0.5        # egg-egg edge probability (q > p)
theta <- runif(n, 0.8, 1.2)
seed <- 123

A_static <- generate_kidney_egg(n, m, p, q, theta = theta, seed = seed)
g_static <- graph_from_adjacency_matrix(A_static, mode = "undirected")

# Color nodes by community
egg_nodes <- which(rowSums(A_static) > quantile(rowSums(A_static), 0.75))
V(g_static)$color <- ifelse(1:n %in% egg_nodes, "tomato", "skyblue")

plot(g_static,
     main = "Static Kidney-Egg Graph",
     vertex.label = NA,
     vertex.size = 20,
     edge.color = "gray50",
     layout = layout_with_fr)

# ----- Test: Dynamic Kidney-Egg Graph -----
cat("Generating dynamic Kidney-Egg graph sequence...\n")

T <- 10
anomaly_start <- 4
anomaly_end <- 6

A_list <- generate_dynamic_kidney_egg(n, m, p, q, T,
                                      anomaly_start, anomaly_end,
                                      theta = theta, seed = seed)

# Plot graphs before, during, and after anomaly
par(mfrow = c(1, 3))  # 3 plots side by side

plot(graph_from_adjacency_matrix(A_list[[3]]),
     main = "Before Anomaly (t = 3)",
     vertex.label = NA,
     vertex.size = 10,
     edge.color = "gray70",
     layout = layout_with_fr)

plot(graph_from_adjacency_matrix(A_list[[5]]),
     main = "During Anomaly (t = 5)",
     vertex.label = NA,
     vertex.size = 20,
     edge.color = "gray50",
     layout = layout_with_fr)

plot(graph_from_adjacency_matrix(A_list[[8]]),
     main = "After Anomaly (t = 8)",
     vertex.label = NA,
     vertex.size = 20,
     edge.color = "gray80",
     layout = layout_with_fr)

# Reset layout
par(mfrow = c(1, 1))
