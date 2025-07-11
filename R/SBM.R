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

# Generates a dynamic SB
#  optional: generate a dynamic DCSBM
#' @export
generate_dynamic_sbm <- function(n, K, Z, B, new_B,theta=NULL, T = 10,
                                 persistence=0.1, start_time,
                                 end_time=start_time,
                                 theta_fluctuate=TRUE,
                                 theta_spread_change = NULL,
                                 theta_spread_blocks = NULL,
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
  # seed: for reproducibility

  # List of adjacency matrices at each time step
  A_list <- list()
  # Starting persistence should be 0 until changepoint occurs
  cur_persistence = 0
  B_t <- B #active B

  if (!is.null(theta) && !is.null(theta_spread_change)) {
    base_theta_min <-min(theta)
    base_theta_max <- max(theta)
  }

  for (t in 1:T) {
    # Changepoint occurs, persistence changes
    if (t == start_time){
      cur_persistence = persistence
      B_t <- new_B
    }
    for (i in 1:n) {
      # Update block membership
      if (runif(1) < persistence) {
        Z[i, ] <- rep(0, K)
        Z[i, sample(1:K, 1)] <- 1
      }
      # Update Degree parameters (random fluctations)
      if (!is.null(theta_fluctuate) && theta_fluctuate) {
        theta <- theta * (1 + runif(n, -0.1, 0.1))
      }

      if (!is.null(theta) &&
          !is.null(theta_spread_change) &&
          !is.null(theta_spread_blocks)) {
        if (t >= start_time && t <= end_time) {
          node_labels <- apply(Z, 1, function(row) which(row == 1))
          for (r in theta_spread_blocks) {
            ids <- which(node_labels == r)
            theta[ids] <- runif(length(ids),
                                min = base_theta_min - theta_spread_change,
                                max = base_theta_max + theta_spread_change)
          }
        }
      }

    }
    A <- generate_sbm(n, K, Z, B_t, theta)
    A_list[[t]] <- A
    # Changepoint ends, persistence back to 0
    if (t == end_time) {
      cur_persistence = 0
      B_t <- B
    }
  }

  return(A_list)
}
