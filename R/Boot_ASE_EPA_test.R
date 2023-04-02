library(Matrix)
library(irlba)


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

# Example: LowRankTests(A, B, r, sig, bs)

#' LowRankTests
#'
#' This function is used to perform the two-sample test to reveal whether the two set of networks from identical population when m=1
#' This function performs the Boot-ASE and Boot-EPA tests, which are two-sample testing procedures for large graphs.
#' This test includes tests which are based on adjacency spectral embedding (ASE) and estimated population adjacency (EPA), respectively. \cr\cr
#' The test involves computing the rejection decision and the p-value, which is computed to ensure that the null hypothesis (i.e., both populations have the same population adjacency) is rejected only when the test statistic is in the upper α-quantile for both approximate distributions.
#'
#'
#' @param A cell array containing network in 1st population; each cell is a sparse EPAacency matrix
#' @param B cell array containing network in 2nd population, each cell being a sparse EPAacency matrix
#' @param r scalar specifying rank of population EPAacency
#' @param sig significance level for acceptance of null hypothesis
#' @param bs bootstrap times
#'
#' @return
#' tASE: output of Boot-ASE test (1 if null is rejected, 0 otherwise)\cr
#' tEPA: output of Boot-EPA test (1 if null is rejected, 0 otherwise)\cr
#' pvalASE: p-value for Boot-ASE test\cr
#' pvalEPA: p-value for Boot-EPA test
#'
#' @export
#'
#' @examples
#' A <- genSparseGraph(1, model1) # sample size, m=1
#' B <- genSparseGraph(1, model2)
#' test_result = LowRankTests(A, B, 100, 0.05, 200)
#'
#'
LowRankTests <- function(A, B, r, sig, bs) {

  # compute the test statistic
  testStat <- Boot_ASE_computeStat(A[[1]], B[[1]], r)
  # bootstrap
  bsStat <- matrix(0, nrow = bs, ncol = 2)

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
    bsStat[b,] <- Boot_ASE_computeStat(C[[1]], C[[2]], r)
  }

  # compute the p-value
  px <- (colSums(bsStat > matrix(testStat, nrow = bs, ncol = 2, byrow = TRUE)) + 0.5) / bs


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
    bsStat[b,] <- Boot_ASE_computeStat(C[[1]], C[[2]], r)
  }

  # compute the p-value
  py <- (colSums(bsStat > matrix(testStat, nrow = bs, ncol = 2, byrow = TRUE)) + 0.5) / bs
  # p-values and acceptance/rejection
  pvalASE <- max(px[1], py[1])
  tASE <- (pvalASE <= sig)

  pvalEPA <- max(px[2], py[2])
  tEPA <- (pvalEPA <= sig)
  return(list(tASE = tASE, tEPA = tEPA, pvalASE = pvalASE, pvalEPA = pvalEPA))
}




# Auxiliary functions


# Boot-ASE/EPA Compute Statistics -----------------------------------------

# compute the test statistic in (7) and (8)
# ASE statistic and EPA statistic


#' Boot-ASE/EPA Compute Statistics
#'
#'
#'
#' This function computes the test statistics for Boot-ASE/EPA test, serving as an auxiliary functions for the function: LowRankTests
#'
#' @param A1 a sparse EPAacency matrix
#' @param B1 a sparse EPAacency matrix
#' @param r scalar specifying rank of population EPAacency
#'
#' @return ASE statistic and EPA statistic
#' @export
#'
#' @examples
#' simu = twosamp_twoblock(n = 4, m = 4, p = 0.3, q = 0.3, epsilon = 0.7)
#' Boot_ASE_computeStat(simu$A_G[[1]], simu$A_H[[1]], 100)
#'
#'
Boot_ASE_computeStat <- function(A1, B1, r) {
  # svd of A1
  u1svd <- svd(as.matrix(A1))
  u1 <- u1svd$u[,1:r]
  s1 <- diag(u1svd$d[1:r])
  v1 <- u1svd$v[,1:r]

  # svd of B1
  u2svd <- svd(as.matrix(B1))
  u2 <- u2svd$u[,1:r]
  s2 <- diag(u2svd$d[1:r])
  v2 <- u2svd$v[,1:r]

  X <- u1 %*% sqrt(s1)
  Y <- u2 %*% sqrt(s2)

  # compute W
  uvSVD <- svd(t(X) %*% Y)
  W <- uvSVD$u %*% uvSVD$v
  stats <- numeric(2)

  # ASE statistic
  # min_W||X-YW||_F solved using orthogonal Procrustes problem, whose solution
  # is Wopt = UV' where X'Y = USV'
  stats[1] <- norm(X - Y %*% W, type = "F")

  # EPA statistic (difference of estimated population EPAacencies)
  stats[2] <- norm(u1 %*% s1 %*% t(v1) - u2 %*% s2 %*% t(v2), type = "F")
  return(stats)
}



