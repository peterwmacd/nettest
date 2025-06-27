#' @export

# Generates a static Kidney-Egg model as a special case of SBM
generate_kidney_egg <- function(n, m, p, q, theta = NULL, seed=NULL) {
  # n: Total # of nodes
  # m: # of "Egg" nodes (chatter nodes)
  # p: Probability of K-K edges and K-E edges
  # q: Probability of E-E edges (should be greater than p, hence the "chatter")
  # theta: (Optional) Degree correction vector of length n
  # seed: (Optional) for reproducibility

  # Output: Adjacency matrix (nxn) of the generated undirected graph

  if (!is.null(seed)) set.seed(seed)

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

#' @export
# Generates a controlled dynamic Kidney-Egg model
generate_dynamic_kidney_egg <- function(n, m, p, q, T = 10,
                                        anomaly_start=5, anomaly_end=7,
                                        theta = NULL, seed=NULL) {
  # # n: Total # of nodes
  # m: # of "Egg" nodes (chatter nodes)
  # p: Probability of K-K edges and K-E edges
  # q: Probability of E-E edges (should be greater than p, hence the "chatter"
  # T: Total # of time steps
  # start: timestep when egg apears
  # end: timestep when egg disappears
  # theta: (Optional) Degree correction vector of length n
  # seed: (Optional) for reproducibility

  # Output: List of Adjacency matrices (nxn) of the generated undirected graphs

  if (!is.null(seed)) set.seed(seed)

  # List of adjecency matrices for each time stepr
  A_list <- list()

  for (t in 1:T) {
    step_seed <- if (!is.null(seed)) seed + t else NULL

    if (t >= anomaly_start && t <= anomaly_end) {
      # Generate Kidney Egg Model during chatter event (anomaly)
      A <- generate_kidney_egg(n, m, p, q, theta, seed=step_seed)
    } else {
      # There are no egg nodes (normal behaviour): regular ER(n,p)
      # Note: ER model has 0 egg nodes and p=q
      A <- generate_kidney_egg(n, 0, p, p, theta, seed=step_seed)
    }

    A_list[[t]] <- A
  }

  return(A_list)
}
