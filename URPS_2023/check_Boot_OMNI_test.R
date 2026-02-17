# Testing OMNI
devtools::load_all()

r <- 4
sig <- 0.05
bs <- 2
num_matrices <- 5

# Generate random adjacency matrices
generate_random_adjacency <- function(n) {
  # Returns a random adjacency matrix of size n x n
  return(matrix(sample(0:10, n * n, replace = TRUE), ncol = n))
}

# Create cell arrays containing multiple adjacency matrices for the first population
A <- lapply(1:num_matrices, function(x) generate_random_adjacency(n = 10))

# Create cell arrays containing multiple adjacency matrices for the second population
B <- lapply(1:num_matrices, function(x) generate_random_adjacency(n = 10))

# Perform OmniRankTests
OmniRankTests(A, B, r, sig, bs)
