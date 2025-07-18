#' Plot Simulation Summary with Shewhart Control Limits and Auto-Scaled Axes
#'
#' @param dynamic_networks List of adjacency matrices (one per time step)
#' @param node_labels Vector of length n indicating community labels (1 or 2)
#' @param sim_title Title to display above the plot
#' @param baseline_window Vector of time steps to use for control limit estimation (default: 1:25)
#'
#' @return A 2x2 base R plot panel with Shewhart-style bounds and visible control limits
#' @export
plot_simulation_summary <- function(dynamic_networks, node_labels, sim_title = "Simulation Summary", baseline_window = 1:25) {

  estimate_block_probs <- function(adj, labels) {
    ids_1 <- which(labels == 1)
    ids_2 <- which(labels == 2)

    A_11 <- adj[ids_1, ids_1]
    A_22 <- adj[ids_2, ids_2]
    A_12 <- adj[ids_1, ids_2]

    p11 <- sum(A_11[lower.tri(A_11)]) / choose(length(ids_1), 2)
    p22 <- sum(A_22[lower.tri(A_22)]) / choose(length(ids_2), 2)
    p12 <- sum(A_12) / (length(ids_1) * length(ids_2))

    return(c(p11 = p11, p12 = p12, p22 = p22))
  }

  compute_control_limits <- function(series, window, k = 3) {
    baseline <- series[window]
    m <- mean(baseline)
    s <- sd(baseline)
    if (is.na(s) || s == 0) {
      return(c(lower = NA, center = m, upper = NA))
    } else {
      return(c(lower = m - k * s, center = m, upper = m + k * s))
    }
  }

  draw_shewhart_plot <- function(series, limits, main_label) {
    y_range <- range(c(series, limits["lower"], limits["upper"]), na.rm = TRUE)
    plot(series, type = "l", ylab = "Value", xlab = "Time", main = main_label, ylim = y_range)
    abline(h = limits["center"], col = "blue")
    if (!is.na(limits["lower"])) abline(h = limits["lower"], col = "red", lty = 2)
    if (!is.na(limits["upper"])) abline(h = limits["upper"], col = "red", lty = 2)
  }

  # Estimate connection probabilities
  block_prob_ts <- t(sapply(dynamic_networks, estimate_block_probs, labels = node_labels))

  # Estimate s1 using Louvain clustering
  pi_max_ts <- sapply(dynamic_networks, function(adj) {
    g <- graph_from_adjacency_matrix(adj, mode = "undirected")
    cl <- cluster_leading_eigen(g)
    sizes <- sizes(cl)
    largest <- which.max(sizes)
    mean(membership(cl) == largest)
  })

  # Extract individual series
  p11_ts <- block_prob_ts[, "p11"]
  p12_ts <- block_prob_ts[, "p12"]
  p22_ts <- block_prob_ts[, "p22"]

  # Compute control limits
  pi_max_lim  <- compute_control_limits(pi_max_ts, baseline_window)
  p11_lim <- compute_control_limits(p11_ts, baseline_window)
  p12_lim <- compute_control_limits(p12_ts, baseline_window)
  p22_lim <- compute_control_limits(p22_ts, baseline_window)

  # Plot 2x2
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  draw_shewhart_plot(pi_max_ts,  pi_max_lim,  expression(pi[max]))
  draw_shewhart_plot(p11_ts, p11_lim, expression(hat(P)[11]))
  draw_shewhart_plot(p12_ts, p12_lim, expression(hat(P)[12]))
  draw_shewhart_plot(p22_ts, p22_lim, expression(hat(P)[22]))
  mtext(sim_title, side = 3, outer = TRUE, line = -1.5, cex = 1.4, font = 2)
  par(mfrow = c(1, 1))
}



