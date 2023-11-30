# Update dependencies
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!requireNamespace("irlba", quietly = TRUE)) {
  install.packages("irlba")
}


genIERGraph <- function(m, model) {
  A <- list()
  for (i in 1:m) {
    EA <- upper.tri(model$P)
    A1 <- matrix(runif(length(EA)) < EA, nrow = ncol(model$P))
    A1 <- A1 + t(A1)  # Ensure symmetry by adding the transpose
    A[[i]] <- A1
  }
  return(A)
}

get_omnibus_matrix_sparse <- function(matrices) {
  rows <- list()
  
  # Iterate over each column
  for (column_index in seq_along(matrices)) {
    current_row <- list()
    
    for (row_index in seq_along(matrices)) {
      if (row_index == column_index) {
        # we are on the diagonal, we do not need to perform any calculation and instead add the current matrix
        # to the current_row
        current_row[[row_index]] <- matrices[[column_index]]
      } else {
        # otherwise we are not on the diagonal and we average the current_matrix with the matrix at row_index
        # and add that to our current_row
        matrices_averaged <- (matrices[[column_index]] + matrices[[row_index]]) * 0.5
        current_row[[row_index]] <- matrices_averaged
      }
    }
    
    result <- matrix(0, nrow=0, ncol=ncol(current_row[[1]]))
    
    # Iterate through the current_row list and vertically concatenate matrices into the result matrix
    for (i in 1:length(current_row)) {
      result <- rbind(result, current_row[[i]])
    }
    
    # row
    rows[[column_index]] <- result
  }
  
  # Combine rows to create the omnibus matrix
  omnibus_matrix <- matrix(0, nrow = nrow(rows[[1]]), ncol=0)
  
  # Iterate through the rows list and horizontally concatenate matrices into the result matrix
  for (i in 1:length(rows)) {
    omnibus_matrix <- cbind(omnibus_matrix, rows[[i]])
  }
  
  return(omnibus_matrix)
}


# OMNI statistic and EPA statistic
Boot_OMNI_computeStat <- function(C, r) {
  # Matrix M
  M <- get_omnibus_matrix_sparse(C)
  
  # Obtain omni
  u1svd <- svd(as.matrix(M))
  u1 <- u1svd$u[, 1:r]
  s1 <- diag(u1svd$d[1:r])
  v1 <- u1svd$v[, 1:r]
  M_omni <- u1 %*% sqrt(s1)
  
  # Select first N rows and last N rows
  n <- nrow(C[[1]])
  A <- M_omni[1:n, , drop = FALSE]
  B <- M_omni[(nrow(M_omni) - n + 1):nrow(M_omni), , drop = FALSE]
  
  # ASE statistic
  # min_W||X-Y||_F solved using orthogonal Procrustes problem, whose solution
  stats <- norm(A - B, type = "F")

  return(stats)
}

# Boot-ASE/EPA Test -------------------------------------------------------

  # Returns acceptance/rejection for bootstrapped tests Boot-ASE (arXiv:1403.7249)
  # Boot-EPA (statistic from arXiv:1606.02401)
  # Note:
  # All graphs are assumed to unweighted, undirected, and defined on a common vertex set.
  # Only the first network from each population is tested
  
  # Input:
  #   A: cell array containing network in 1st population; each cell is a sparse EPAacency matrix
  #   B: cell array containing network in 2nd population, each cell being a sparse EPAacency matrix
  #   r = scalar specifying rank of population EPAacency
  #   sig: significance level for acceptance of null hypothesis
  #   bs: bootstrap times
  
  # Output:
  #   tASE: output of Boot-ASE test (1 if null is rejected, 0 otherwise)
  #   tEPA: output of Boot-EPA test (1 if null is rejected, 0 otherwise)
  #   pvalASE: p-value for Boot-ASE test
  #   pvalEPA: p-value for Boot-EPA test
  
  # Example: OmniRankTests(A, B, r, sig, bs)

OmniRankTests <- function(A, B, r, sig, bs) {
  # compute the test statistic
  # svd of A1
  u1svd <- svd(as.matrix(A[[1]]))
  u1 <- u1svd$u[,1:r]
  s1 <- diag(u1svd$d[1:r])
  v1 <- u1svd$v[,1:r]
  
  # svd of B1
  u2svd <- svd(as.matrix(B[[1]]))
  u2 <- u2svd$u[,1:r]
  s2 <- diag(u2svd$d[1:r])
  v2 <- u2svd$v[,1:r]
  
  X <- u1 %*% sqrt(s1)
  Y <- u2 %*% sqrt(s2)
  
  # compute 
  testStat <- norm(X - Y, type = "F")
  
  # bootstrap
  bsStat <- matrix(0, nrow = bs, ncol = 1)
  
  # generate samples from E[A]
  usv <- svd(as.matrix(A[[1]]), r, r)
  u <- usv$u
  s <- diag(usv$d[1:r])
  v <- usv$v
  EA <- u %*% s %*% t(v)
  EA <- EA - diag(diag(EA))
  EA1 <- ifelse(EA > 0, EA, 0)
  EA <- ifelse(EA1 > 1, 1, EA1)
  model <- list(name = "IER", n = nrow(EA), P = EA)
  
  for (b in 1:bs) {
    C <- genIERGraph(2, model)
    bsStat[b,] <- Boot_OMNI_computeStat(C, r)
  }
  
  # compute the p-value
  px <- (colSums(bsStat > matrix(testStat, nrow = bs, ncol = 1, byrow = TRUE)) + 0.5) / bs
  
  # generate samples from E[B]
  usv <- svd(as.matrix(B[[1]]), r, r)
  u <- usv$u
  s <- diag(usv$d[1:r])
  v <- usv$v
  EB <- u %*% s %*% t(v)
  EB <- EB - diag(diag(EB))
  EB1 <- ifelse(EB > 0, EB, 0)
  EB <- ifelse(EB1 > 1, 1, EB1)
  model$P <- EB
  
  for (b in 1:bs) {
    C <- genIERGraph(2, model)
    bsStat[b,] <- Boot_OMNI_computeStat(C, r)
  }
  
  # compute the p-value
  py <- (colSums(bsStat > matrix(testStat, nrow = bs, ncol = 1, byrow = TRUE)) + 0.5) / bs
  # p-values and acceptance/rejection
  pvalOMNI <- max(px[1], py[1])
  tOMNI <- ( pvalOMNI <= sig)
  return(list(tOMNI = tOMNI, pvalOMNI = pvalOMNI))
}

# Testing
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







