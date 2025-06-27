devtools::load_all()


# This function generates a static Stochastic Block Model using Datatables
#.  to make an edgelist
#   Note: better for larger datasets
#  optional: generates a Degree Corrected SBM
generate_sbm_dt <- function(n, K, Z, B, theta=NULL) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K)
  # B: Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N)

  # Output: An undirected datatable

  edgelist <- data.table(i=integer(0), j=integer(0))


  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      i_block <- which(Z[i, ] == 1)
      j_block <- which(Z[j, ] == 1)

      prob_ij <- B[i_block, j_block]

      # Degree Correction
      if (!is.null(theta)) {
        prob_ij <- prob_ij * theta[i] * theta[j]
      }

      # Cap probabilities at 1 (in case of degree correction above 1)
      prob_ij <- min(prob_ij, 1)

      if (runif(1) < prob_ij) {
        edgelist <- rbind(edgelist, data.table(i=i,j=j))
      }
    }
  }
  return(edgelist)
}

# Generates a dynamic SB using a datatable
#  optional: generate a dynamic DCSBM
generate_dynamic_sbm_dt <- function(n, K, Z, B, theta=NULL, T = 10,
                                 persistence=0.1) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K)
  # B: Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N)
  # T: # of time steps for dynamic simulation
  # persistence: probability of block change for a node

  # List of adjacency matrices at each time step
  edgelist <- data.table(from=integer(0), to=integer(0), times=vector("list", 0))

  for (t in 1:T) {
    for (i in 1:n) {
      if (runif(1) < persistence) {
        # Randomly reassign to a new block
        new_block <- sample(1:K, 1)
        Z[i, ] <- 0
        Z[i, new_block] <- 1
      }
    }

    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        i_block <- which(Z[i, ] == 1)
        j_block <- which(Z[j, ] == 1)

        prob_ij <- B[i_block, j_block]

        # Degree Correction
        if (!is.null(theta)) {
          prob_ij <- prob_ij * theta[i] * theta[j]
        }

        # Cap probabilities at 1 (in case of degree correction above 1)
        prob_ij <- min(prob_ij, 1)

        if (runif(1) < prob_ij) {
          # check if edge exists already
          idx <- edgelist[from==i & to==j, which=TRUE]
          if (length(idx) == 0) {
            # new edge
            edgelist <-rbind(edgelist, data.table(from=i,to=j,
                                                  times=list(c(t))))
          } else {
            # Existing edge, append time step
            edgelist[idx, times := list(c(times[[1]], t))]
          }
        }
      }
    }
  }

  return(edgelist)
}
