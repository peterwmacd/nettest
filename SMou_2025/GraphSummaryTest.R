# Simulate a toy dynamic DCSBM
set.seed(123)
n <- 20
K <- 2
Z <- matrix(0, n, K)
Z[1:(n/2), 1] <- 1
Z[((n/2)+1):n, 2] <- 1

B <- matrix(c(0.2, 0.05, 0.05, 0.2), 2, 2)
theta <- runif(n, 0.8, 1.2)
T <- 30

# Simulate dynamic network
dynamic_networks <- generate_dynamic_sbm(
  n = n, K = K, Z = Z,
  B = B, new_B = B,
  theta = theta,
  T = T,
  persistence = 0,
  start_time = 15,
  end_time = 15
)

# Compute F1–F9 summary statistics over time
F_time_series <- compute_all_F_series(dynamic_networks)

# Plot all nine time series in a 3x3 grid
par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))  # 3x3 layout with tighter margins
for (i in 1:9) {
  plot(F_time_series[[i]], type = "l",
       main = paste0("F", i),
       xlab = "Time", ylab = paste0("F", i))
}
par(mfrow = c(1, 1))
