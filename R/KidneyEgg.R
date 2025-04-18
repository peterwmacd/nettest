library(igraph)
library(Matrix)
source("SBM.R")

# Generates a static Kidney-Egg model as a special case of SBM
generate_kidney_egg <- function(n, m, p, q, theta = NULL) {
  # n: Total # of nodes
  # m: # of "Egg" nodes (chatter nodes)
  # p: Probability of K-K edges and K-E edges
  # q: Probability of E-E edges (should be greater than p, hence the "chatter")
  # theta: (Optional) Degree correction vector of length n
  
  # Output: Adjacency matrix (nxn) of the generated undirected graph
  
  
  
  # Number of blocks: 1= Kidney, 2 = Egg
  K <- 2
  
  # Assign nodes to blocks
  egg_nodes <- sample(1:n, m)
  kidney_nodes <- setdiff(1:n, egg_nodes)
  
  # Creates block membership matrix Z (nxK)
  Z <- matrix(0, nrow = n, ncol = K)
  Z[kidney_nodes, 1] <- 1 # Assign Kidney nodes to block 1
  Z[egg_nodes, 2] <- 1 # Assign Egg nodes to block 2
  
  # Define block connection matrix: B (KxK)
  B <- matrix(p, nrow = K, ncol=K)
  B[2, 2] <- q
  
  A <- generate_sbm(n, K, Z, B, theta)
  
  return(A)
}