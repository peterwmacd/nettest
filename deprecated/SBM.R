
#' Generate a Static Stochastic Block Model (SBM)
#'
#' @description
#' Generates an adjacency matrix from a stochastic block model (SBM).
#' Optionally, degree-corrected SBM (DCSBM) can be simulated by supplying
#' node-specific degree parameters.
#'
#' @param n Integer. Number of nodes in the network.
#' @param K Integer. Number of communities (blocks).
#' @param Z Numeric matrix of size \eqn{n \times K}. Block membership matrix,
#'   where each row is a one-hot vector indicating the community assignment
#'   of a node.
#' @param B Numeric \eqn{K \times K} matrix. Block connection probabilities
#'   between communities.
#' @param theta Optional numeric vector of length \eqn{n}. Degree correction
#'   parameters. If supplied, generates a DCSBM. If omitted, a standard SBM
#'   is generated.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A binary adjacency matrix of size \eqn{n \times n}.
#'
#' @examples
#' # Example: 6 nodes, 2 communities
#' n <- 6; K <- 2
#' Z <- matrix(0, n, K)
#' Z[1:3, 1] <- 1
#' Z[4:6, 2] <- 1
#'
#' B <- matrix(c(0.8, 0.2,
#'               0.2, 0.6), nrow = K, byrow = TRUE)
#'
#' # Generate SBM
#' A <- generate_sbm(n, K, Z, B, seed = 42)
#'
#' # Degree-corrected SBM
#' theta <- runif(n, 0.8, 1.2)
#' A_dc <- generate_sbm(n, K, Z, B, theta = theta, seed = 42)
#'
#' @export
generate_sbm <- function(n, K, Z, B, theta = NULL, seed = NULL) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K)
  # B: Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N)
  # seed: (Optional) Seed for reproducibility

  if (!is.null(seed)) set.seed(seed)

  # Create adjacency matrix
  A <- matrix(0, nrow=n, ncol=n)

  # Determine Edge probabilities
  for (i in 1:n) {
    for (j in i:n) {
      # Get block memberships of nodes i and j
      i_block <- which(Z[i, ] == 1)
      j_block <- which(Z[j, ] == 1)

      # Calculating edge probability using B
      prob_ij <- B[i_block, j_block]
      # Degree Correction (if given theta)
      if (!is.null(theta)) {
        prob_ij <- prob_ij * theta[i] * theta[j]
        # Prevents overflow
        prob_ij <- min(prob_ij, 1)
      }

      # Determine (in)existence of edge
      if (i != j && runif(1) < prob_ij) {
        A[i,j] <- 1
        A[j,i] <- 1
      }
    }
  }

  return(A)
}

#' Generate a Dynamic Stochastic Block Model (SBM)
#'
#' @description
#' Simulates a sequence of adjacency matrices over time under a dynamic SBM or
#' dynamic degree-corrected SBM (DCSBM). Supports structural changes such as
#' changepoints, community merging, and community splitting.
#'
#' @param n Integer. Number of nodes.
#' @param K Integer. Number of communities at initialization.
#' @param Z Numeric matrix of size \eqn{n \times K}. Initial block membership
#'   matrix (rows are one-hot vectors).
#' @param B Numeric \eqn{K \times K} matrix. Baseline block connection probabilities.
#' @param new_B Numeric \eqn{K \times K} matrix. Block connection probabilities
#'   to use during changepoint periods.
#' @param theta Optional numeric vector of length \eqn{n}. Degree correction
#'   parameters. If supplied, generates a DCSBM.
#' @param T Integer. Number of time steps to simulate.
#' @param persistence Numeric in [0,1]. Probability that a node changes its
#'   community label at each time step during a changepoint.
#' @param start_time Integer. First time step of the changepoint.
#' @param end_time Integer. Last time step of the changepoint. Defaults to
#'   \code{start_time}.
#' @param theta_fluctuate Logical. If TRUE, node degree parameters fluctuate
#'   randomly over time.
#' @param theta_spread_change Optional numeric. If provided, widens the range
#'   of \code{theta} for specified communities during a changepoint.
#' @param theta_spread_blocks Optional integer vector. Community indices to
#'   which \code{theta_spread_change} is applied.
#' @param merge_communities Logical. If TRUE, all communities are merged into
#'   a single community at \code{start_time}.
#' @param split_community Logical. If TRUE, community 1 is split into two
#'   subcommunities at \code{start_time}.
#' @param split_within Numeric. Probability of edges within new subcommunities
#'   after a split.
#' @param split_between Numeric. Probability of edges between new subcommunities
#'   after a split.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{adj_list} – List of adjacency matrices (length \code{T}).
#'   \item \code{Z_list} – List of membership matrices (length \code{T}).
#' }
#'
#' @examples
#' # Example: dynamic SBM with local changepoint
#' n <- 10; K <- 2; T <- 20
#' Z <- matrix(0, n, K)
#' Z[1:(n/2), 1] <- 1
#' Z[(n/2+1):n, 2] <- 1
#'
#' B <- matrix(c(0.3, 0.1,
#'               0.1, 0.3), nrow = K, byrow = TRUE)
#' new_B <- B; new_B[1,1] <- 0.5
#'
#' sim <- generate_dynamic_sbm(
#'   n = n, K = K, Z = Z,
#'   B = B, new_B = new_B,
#'   T = T, start_time = 10, end_time = 20,
#'   persistence = 0, theta_fluctuate = FALSE
#' )
#'
#' length(sim$adj_list)  # 20 adjacency matrices
#'
#' @export

generate_dynamic_sbm <- function(n, K, Z, B, new_B,theta=NULL, T = 10,
                                 persistence=0.1, start_time,
                                 end_time=start_time,
                                 theta_fluctuate=TRUE,
                                 theta_spread_change = NULL,
                                 theta_spread_blocks = NULL,
                                 merge_communities = FALSE,
                                 split_community = FALSE,
                                 split_within = 0.2,
                                 split_between = 0.15,
                                 seed = NULL) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K)
  # B: Block connection matrix (K x K)
  # new_B: New Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N)
  # T: # of time steps for dynamic simulation
  # persistence: probability of block change for a node
  # start_time: starting timestep of the changepoint
  # end_time: last timestep of the changepoint
  # theta_fluctuate: determines if theta should change over time
  # theta_spread_change: Optional value to widen theta range during changepoint
  # theta_spread_blocks: Optional vector of community indices to apply spread change to
  # merge_communities: determines whether to merge communities
  # split_community: determines whether to split a community
  # split_within: controls within-community probability
  # split_between: controls between-community probability
  # seed: for reproducibility

  # List of adjacency matrices at each time step
  A_list <- list()
  # List of membership matrices over time
  Z_list <- list()
  # Starting persistence should be 0 until changepoint occurs
  cur_persistence = 0
  B_t <- B #active B
  cur_K <- K

  if (!is.null(theta) && !is.null(theta_spread_change)) {
    base_theta_min <-min(theta)
    base_theta_max <- max(theta)
  }

  for (t in 1:T) {
    for (i in 1:n) {
      # Update block membership
      if (runif(1) < persistence) {
        Z[i, ] <- rep(0, K)
        Z[i, sample(1:K, 1)] <- 1
      }
    }

    # Update Degree parameters (random fluctations)
    if (isTRUE(theta_fluctuate)) {
      theta <- theta * (1 + runif(n, -0.1, 0.1))
    }

    if (!is.null(theta) &&
        !is.null(theta_spread_change) &&
        !is.null(theta_spread_blocks) &&
        t >= start_time && t <= end_time) {
      node_labels <- apply(Z, 1, function(row) which(row == 1))
      for (r in theta_spread_blocks) {
        ids <- which(node_labels == r)
        theta[ids] <- runif(length(ids),
                            min = base_theta_min - theta_spread_change,
                            max = base_theta_max + theta_spread_change)
      }
    }

    if (t == start_time) {
      # Changepoint occurs, persistence changes
      cur_persistence = persistence
      B_t <- new_B

      # Merge?
      if (merge_communities) {
        K <- 1
        Z <- matrix(1, nrow=n, ncol = 1)
        B_t <- matrix(mean(B_t), 1, 1)
      }
      # Split?
      if (split_community) {
        idx_comm1 <- which(apply(Z, 1, function(row) which(row == 1)) == 1)
        split_A <- idx_comm1[1:floor(length(idx_comm1)/2)]
        split_B <- setdiff(idx_comm1, split_A)

        Z_new <- matrix(0, n, K + 1)
        Z_new[, 2:(K+1)] <- Z  # shift existing groups to 2 and 3
        Z_new[split_A, 1] <- 1
        Z_new[split_B, 2] <- 1

        Z <- Z_new
        K <- K + 1

        # Expand B matrix
        B_t <- matrix(split_between, K, K)
        diag(B_t) <- split_within
      }
    }


    A <- generate_sbm(n, K, Z, B_t, theta)
    A_list[[t]] <- A
    Z_list[[t]] <- matrix(Z, nrow = n)

    # Changepoint ends, persistence back to 0
    if (t == end_time) {
      cur_persistence = 0
      B_t <- B
    }
  }

  return(list(adj_list = A_list, Z_list = Z_list))
}
