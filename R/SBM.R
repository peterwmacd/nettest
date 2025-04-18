library(igraph)
library(Matrix)


# This function generates a static Stochastic Block Model
#  optional: generates a Degree Corrected SBM
generate_sbm <- function(n, K, Z, B, theta = NULL) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K) 
  # B: Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N) 
  
  
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
generate_dynamic_sbm <- function(n, K, Z, B, theta=NULL, T = 10,
                                 persistence=0.1) {
  # n: # of nodes
  # K: # of blocks (communities)
  # Z: Block membership matrix (n x K) 
  # B: Block connection matrix (K x K)
  # theta: (Optional, if DCSBM is desired) Vector of degree parameters (length N) 
  # T: # of time steps for dynamic simulation
  # persistence: probability of block change for a node
  
  # List of adjacency matrices at each time step
  A_list <- list()
  
  for (t in 1:T) {
    for (i in 1:n) {
      # Update block membership
      if (run_if(1) < persistence) {
        Z[i, ] <- 0
        Z[i, sample(1:K, 1)] <- 1
      }
      # Update Degree parameters (random fluctations)
      if (!is.null(theta)) {
        theta <- theta * (1 + runif(n, -0.1, 0.1))
      }
    }
    A <- generate_sbm(n, K, Z, B, theta)
    A_list[t] = A
  }
  
  return(A_list)
}