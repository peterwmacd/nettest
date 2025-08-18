#' @export

# This function generates a static Stochastic Block Model
#  optional: generates a Degree Corrected SBM
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

# Generates a dynamic SBM
#  optional: generate a dynamic DCSBM
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
    if (!is.null(theta_fluctuate)) {
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
