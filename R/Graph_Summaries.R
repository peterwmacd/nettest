#' Compute F1: Number of edges in the graph
#'
#' @param adj Adjacency matrix (numeric, symmetric, binary)
#'
#' @return Total number of edges (numeric)
#' @export
compute_F1 <- function(adj) {
  sum(adj) / 2
}

#' Compute F2: Maximum degree
#'
#' @param adj Adjacency matrix (numeric, symmetric, binary)
#'
#' @return Maximum node degree (numeric)
#' @export
compute_F2 <- function(adj) {
  max(rowSums(adj))
}


#' Compute F3: Largest eigenvalue of the adjacency matrix
#'
#' @param adj Adjacency matrix (numeric, symmetric)
#'
#' @return Largest eigenvalue (numeric), as proxy for max average degree
#' @export
compute_F3 <- function(adj) {
  max(eigen(adj, only.values = TRUE)$values)
}


#' Compute scan statistic for a given neighborhood scale
#'
#' @param g An igraph graph object
#' @param k Integer, neighborhood order (e.g., 1 for 1-hop)
#'
#' @return Maximum number of edges in any k-hop neighborhood (integer)
#' @export
compute_scan_stat <- function(g, k) {
  max(sapply(V(g), function(v) {
    nbhd <- neighborhood(g, order = k, nodes = v)[[1]]
    subg <- induced_subgraph(g, nbhd)
    ecount(subg)
  }))
}


#' Compute F4-6: Scan statistics (1,2,3-hop neighborhood)
#'
#' @param adj Adjacency matrix
#'
#' @return Maximum number of edges in a neighbourhood of a given size
#' @export
compute_F4 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  compute_scan_stat(g, 1)
}

#' @export
#' @rdname compute_F4
compute_F5 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  compute_scan_stat(g, 2)
}

#' @export
#' @rdname compute_F4
compute_F6 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  compute_scan_stat(g, 3)
}


#' Compute F7: Triangle count
#'
#' @param adj Adjacency matrix
#'
#' @return Total number of triangles in the graph
#' @export
compute_F7 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  sum(count_triangles(g))
}


#' Compute F8: Global clustering coefficient
#'
#' @param adj Adjacency matrix
#'
#' @return Global transitivity (numeric)
#' @export
compute_F8 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  transitivity(g, type = "global")
}


#' Compute F9: Negated average path length (reachability metric)
#'
#' @param adj Adjacency matrix
#'
#' @return Negative of average shortest path length (numeric)
#' @export
compute_F9 <- function(adj) {
  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  if (is_connected(g)) {
    return(-mean_distance(g))
  } else {
    dists <- distances(g)
    n <- vcount(g)
    max_dist <- max(dists[is.finite(dists)])
    dists[!is.finite(dists)] <- 2 * max_dist
    return(-mean(dists[lower.tri(dists)]))
  }
}

# Wrapper to compute all
compute_all_F <- function(adj) {
  list(
    F1 = compute_F1(adj),
    F2 = compute_F2(adj),
    F3 = compute_F3(adj),
    F4 = compute_F4(adj),
    F5 = compute_F5(adj),
    F6 = compute_F6(adj),
    F7 = compute_F7(adj),
    F8 = compute_F8(adj),
    F9 = compute_F9(adj)
  )
}

#' Compute time series of F1–F9 statistics across a dynamic network
#'
#' @param adj_list A list of adjacency matrices representing a dynamic network sequence.
#'
#' @return A data frame with one row per time step and columns F1 through F9.
#' Each entry contains the value of the corresponding summary statistic at that time.
#' @export
compute_all_F_series <- function(adj_list) {
  T <- length(adj_list)

  F_matrix <- matrix(NA, nrow = T, ncol = 9)
  colnames(F_matrix) <- paste0("F", 1:9)

  for (t in 1:T) {
    adj <- adj_list[[t]]
    F_matrix[t, ] <- unlist(compute_all_F(adj))
  }

  return(as.data.frame(F_matrix))
}

#' Plot F1–F9 with Shewhart-style control bands
#'
#' @param F_time_series A data frame with 9 columns named F1 to F9, each a time series
#' @param baseline_window Vector of indices (default 1:25) used to estimate baseline control limits
#' @param sim_title Optional: title above the 3x3 grid
#'
#' @return None. Displays a 3×3 plot of F-statistics with control bands
#' @export
plot_F_summary_with_control_bands <- function(F_time_series, baseline_window = 1:25, sim_title = "") {
  stopifnot(ncol(F_time_series) == 9)
  par(mfrow = c(3, 3), mar = c(4, 4, 2, 1))

  for (j in 1:9) {
    series <- F_time_series[[j]]
    name <- names(F_time_series)[j]

    baseline <- series[baseline_window]
    m <- mean(baseline)
    s <- sd(baseline)

    UCL <- m + 3 * s
    LCL <- m - 3 * s
    y_range <- range(c(series, UCL, LCL), na.rm = TRUE)

    plot(series, type = "l",
         ylim = y_range,
         ylab = "Value", xlab = "Time",
         main = name)
    abline(h = m, col = "blue")
    abline(h = c(UCL, LCL), col = "red", lty = 2)
  }

  mtext(sim_title, side = 3, outer = TRUE, line = -1.5, font = 2, cex = 1.4)
  par(mfrow = c(1, 1))
}


