#' Simulate dynamic SBM with structural events (splits or merges at start_time)
#'
#' @param n Number of nodes.
#' @param K Initial number of communities.
#' @param Z Initial block membership matrix (n × K).
#' @param B Initial block connectivity matrix (K × K).
#' @param theta Optional degree correction vector (length n).
#' @param T Total number of time steps.
#' @param merge Logical; if TRUE, all communities merge into one at start_time.
#' @param split Logical; if TRUE, community 1 is split into two at start_time.
#' @param split_within Within-community edge probability post-split.
#' @param split_between Between-community edge probability post-split.
#' @param persistence Probability a node changes community during changepoint.
#' @param start_time Time step when changepoint + structural event happens.
#' @param end_time Last time step of changepoint.
#' @param theta_fluctuate Logical, whether theta fluctuates over time.
#' @param theta_spread_change Amount to expand theta range during changepoint.
#' @param theta_spread_blocks Indices of communities affected by theta spread.
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with:
#' \describe{
#'   \item{adj_list}{List of adjacency matrices over time.}
#'   \item{Z_list}{List of community membership matrices over time.}
#' }
#' @export
simulate_structural_sbm <- function(n, K, Z, B,new_B=B,
                                    theta = NULL, T=10,
                                    merge = FALSE,
                                    split = FALSE,
                                    split_within = 0.2,
                                    split_between = 0.15,
                                    persistence = 0.1,
                                    start_time,
                                    end_time = start_time,
                                    theta_fluctuate = TRUE,
                                    theta_spread_change = NULL,
                                    theta_spread_blocks = NULL,
                                    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  adj_list <- list()
  Z_list <- list()

  current_Z <- Z
  current_B <- B
  current_K <- K

  for (t in 1:T) {
    # === Merge at changepoint ===
    if (merge && t == start_time) {
      message("Merging communities at t = ", t)
      current_K <- 1
      current_Z <- matrix(1, nrow = n, ncol = 1)
      current_B <- matrix(mean(current_B), 1, 1)
    }

    # === Split at changepoint ===
    if (split && t == start_time) {
      message("Splitting community 1 at t = ", t)
      comm1_ids <- which(max.col(current_Z) == 1)
      split_A <- comm1_ids[1:floor(length(comm1_ids) / 2)]
      split_B <- setdiff(comm1_ids, split_A)

      Z_new <- matrix(0, n, current_K + 1)
      Z_new[, 2:(current_K + 1)] <- current_Z
      Z_new[split_A, 1] <- 1
      Z_new[split_B, 2] <- 1

      current_Z <- Z_new
      current_K <- current_K + 1

      current_B <- matrix(split_between, current_K, current_K)
      diag(current_B) <- split_within
    }

    # === Simulate one time step ===
    single_step <- generate_dynamic_sbm(
      n = n,
      K = current_K,
      Z = current_Z,
      B = current_B,
      theta = theta,
      T = 1,
      persistence = persistence,
      start_time = start_time,
      end_time = end_time,
      new_B = current_B,
      theta_fluctuate = theta_fluctuate,
      theta_spread_change = theta_spread_change,
      theta_spread_blocks = theta_spread_blocks,
      seed = NULL
    )

    adj_list[[t]] <- single_step$adj_list[[1]]
    Z_list[[t]] <- single_step$Z_list[[1]]
    current_Z <- single_step$Z_list[[1]]
  }

  return(list(adj_list = adj_list, Z_list = Z_list))
}
