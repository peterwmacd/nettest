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

#' @export
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
    C <- genSparseGraph(2, model)
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
    C <- genSparseGraph(2, model)
    bsStat[b,] <- Boot_OMNI_computeStat(C, r)
  }

  # compute the p-value
  py <- (colSums(bsStat > matrix(testStat, nrow = bs, ncol = 1, byrow = TRUE)) + 0.5) / bs
  # p-values and acceptance/rejection
  pvalOMNI <- max(px[1], py[1])
  tOMNI <- ( pvalOMNI <= sig)
  return(list(tOMNI = tOMNI, pvalOMNI = pvalOMNI))
}







