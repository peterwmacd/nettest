n <- 50
K <- 2
Z <- matrix(0, n, K)
Z[1:(n/2), 1] <- 1
Z[((n/2)+1):n, 2] <- 1

B <- matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)      # Base connectivity
new_B <- B
new_B[1, 1] <- 0.3                           # Local increase for community 1

theta <- runif(n, 0.8, 1.2)
T <- 50
t_star <- 30
persistence <- 0                            # No label switching


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
  end_time = T,
  theta_fluctuate = FALSE
)

# === Compute F1–F9 summaries ===
F_time_series <- compute_all_F_series(dynamic_networks)

par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))
for (i in 1:9) {
  plot(F_time_series[[i]], type = "l",
       main = paste0("F", i),
       xlab = "Time", ylab = paste0("F", i))
  abline(v = t_star, col = "red", lty = 2)
}

par(mfrow = c(1, 1))
