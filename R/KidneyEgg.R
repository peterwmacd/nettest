#' Generate a Static Kidney–Egg Network (SBM Special Case)
#'
#' @description
#' Builds a “Kidney–Egg” graph as a 2-block Stochastic Block Model (SBM),
#' where block 1 (“Kidney”) is the background and block 2 (“Egg”) is a
#' dense chatter group. Edges within the Egg (2–2) occur with probability \eqn{q},
#' while all other pairs (1–1 and 1–2) occur with probability \eqn{p}.
#' Optionally simulates a degree-corrected SBM (DCSBM) when \code{theta} is supplied.
#'
#' @param n Integer. Total number of nodes.
#' @param m Integer. Number of “Egg” (chatter) nodes; must satisfy \code{0 <= m <= n}.
#' @param p Numeric in [0, 1]. Baseline edge probability for Kidney–Kidney and Kidney–Egg edges.
#' @param q Numeric in [0, 1]. Egg–Egg (block 2–2) edge probability (typically \code{q > p}).
#' @param theta Optional numeric vector of length \code{n}. Degree correction parameters for DCSBM.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A binary \eqn{n \times n} adjacency matrix.
#'
#' @details
#' The model is constructed as a 2-block SBM with membership matrix \eqn{Z} and
#' block matrix \eqn{B}, where \eqn{B_{11} = p}, \eqn{B_{12} = p}, \eqn{B_{22} = q}.
#'
#' @examples
#' set.seed(1)
#' n <- 20; m <- 6
#' A <- generate_kidney_egg(n = n, m = m, p = 0.05, q = 0.25, seed = 123)
#' # igraph::plot(igraph::graph_from_adjacency_matrix(A, mode = "undirected"))
#'
#' # With degree correction
#' theta <- runif(n, 0.8, 1.2)
#' A_dc <- generate_kidney_egg(n, m, p = 0.05, q = 0.25, theta = theta, seed = 123)
#'
#' @seealso \code{\link{generate_sbm}}, \code{\link{generate_dynamic_kidney_egg}}
#' @export
generate_kidney_egg <- function(n, m, p, q, theta = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Basic validation (lightweight to avoid surprises)
  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("n must be a positive integer.")
  if (!is.numeric(m) || length(m) != 1 || m < 0 || m > n) stop("m must be an integer in [0, n].")
  if (!is.numeric(p) || p < 0 || p > 1) stop("p must be in [0, 1].")
  if (!is.numeric(q) || q < 0 || q > 1) stop("q must be in [0, 1].")
  if (!is.null(theta) && length(theta) != n) stop("theta must be length n when supplied.")

  # 2-block SBM: 1 = Kidney, 2 = Egg
  K <- 2

  # Sample Egg and Kidney indices
  egg_nodes <- if (m > 0) sample.int(n, m) else integer(0)
  kidney_nodes <- setdiff(seq_len(n), egg_nodes)

  # Membership matrix Z (n x K)
  Z <- matrix(0, nrow = n, ncol = K)
  if (length(kidney_nodes)) Z[kidney_nodes, 1] <- 1
  if (length(egg_nodes))    Z[egg_nodes, 2]    <- 1

  # Block matrix B (K x K)
  B <- matrix(p, nrow = K, ncol = K)
  B[2, 2] <- q

  # Delegate to SBM generator
  generate_sbm(n = n, K = K, Z = Z, B = B, theta = theta, seed = NULL)
}


#' Generate a Dynamic Kidney–Egg Network (Controlled Window)
#'
#' @description
#' Produces a time series of Kidney–Egg graphs. During an “anomaly window”
#' (\code{anomaly_start} to \code{anomaly_end}), the Egg (size \code{m}) is present
#' with probabilities \code{p} (for 1–1 and 1–2 edges) and \code{q} (for 2–2 edges).
#' Outside the window, the network reverts to an Erdős–Rényi-like background
#' by setting \code{m = 0} and \code{q = p}.
#'
#' @param n Integer. Total number of nodes.
#' @param m Integer. Number of Egg nodes during the anomaly window.
#' @param p Numeric in [0, 1]. Baseline probability for Kidney–Kidney and Kidney–Egg edges.
#' @param q Numeric in [0, 1]. Egg–Egg probability during anomaly window.
#' @param T Integer. Number of time steps (length of the series). Default \code{10}.
#' @param anomaly_start Integer. First time step (inclusive) where the Egg is active.
#' @param anomaly_end Integer. Last time step (inclusive) where the Egg is active.
#' @param theta Optional numeric vector of length \code{n}. Degree correction parameters.
#' @param seed Optional integer. Base random seed for reproducibility.
#'
#' @return A list of \eqn{T} adjacency matrices (each \eqn{n \times n}).
#'
#' @examples
#' set.seed(1)
#' A_list <- generate_dynamic_kidney_egg(
#'   n = 30, m = 8, p = 0.04, q = 0.2,
#'   T = 20, anomaly_start = 8, anomaly_end = 14, seed = 99
#' )
#' length(A_list)   # 20 graphs
#'
#' @seealso \code{\link{generate_kidney_egg}}, \code{\link{generate_sbm}}
#' @export
generate_dynamic_kidney_egg <- function(n, m, p, q, T = 10,
                                        anomaly_start = 5, anomaly_end = 7,
                                        theta = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Basic validation
  if (!is.numeric(T) || length(T) != 1 || T <= 0) stop("T must be a positive integer.")
  if (!is.numeric(anomaly_start) || !is.numeric(anomaly_end)) {
    stop("anomaly_start and anomaly_end must be integers.")
  }
  if (anomaly_start < 1 || anomaly_end < anomaly_start || anomaly_end > T) {
    stop("Require 1 <= anomaly_start <= anomaly_end <= T.")
  }
  if (!is.null(theta) && length(theta) != n) stop("theta must be length n when supplied.")

  A_list <- vector("list", T)

  for (t in seq_len(T)) {
    step_seed <- if (!is.null(seed)) (as.integer(seed) + t) else NULL

    if (t >= anomaly_start && t <= anomaly_end) {
      # Active Egg: m nodes in block 2, q controls 2–2 density
      A_list[[t]] <- generate_kidney_egg(n, m, p, q, theta = theta, seed = step_seed)
    } else {
      # No Egg: set m = 0 and q = p (ER-like background)
      A_list[[t]] <- generate_kidney_egg(n, 0, p, p, theta = theta, seed = step_seed)
    }
  }

  A_list
}
