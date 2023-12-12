# Boot-Spectral Test ------------------------------------------------------

# Permutation based bootstrapped variants of tests in arxiv 1707.00833 (Boot-Frobenius and Boot-Spectral)
# Note:
# All graphs are assumed to unweighted, undirected, and defined on a common vertex set.
# Sample size in each population must be at least 2

# Input:
#   A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#   B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#   sig: significance level for acceptance of null hypothesis
#   bs: bootstrap times

# Output:
#   tFro: output of Boot-Frobenius test (1 if null is rejected, 0 otherwise)
#   tOp: output of Boot-Spectral test (1 if null is rejected, 0 otherwise)
#   pvalFro: p-value for Boot-Frobenius test
#   pvaOP: p-value for Boot-Spectral test

# Example: ShufflingTests(A, B, sig, bs)


#' Shuffling Tests
#'
#'
#' This function performs the Boot-Frobenius and Boot-Spectral tests, which are two-sample testing procedures for large graphs.
#' The Boot-Frobenius test is a two-sample testing procedure for large graphs. The test is based on the Frobenius norm and uses bootstrapping to approximate the null distribution of a test statistic.
#' The Boot-Spectra is based on adjacency spectral embedding (ASE) and uses bootstrapping to approximate the null distribution of a test statistic. \cr\cr
#' The test involves computing the rejection decision and the p-value, which is computed to ensure that the null hypothesis (i.e., both populations have the same population adjacency) is rejected only when the test statistic is in the upper α-quantile for both approximate distributions.
#'
#'
#' @param A cell array containing network in 1st population; each cell is a sparse adjacency matrix
#' @param B cell array containing network in 2nd population, each cell being a sparse adjacency matrix
#' @param sig significance level for acceptance of null hypothesis
#' @param bs bootstrap times
#'
#' @return
#' tFro: output of Boot-Frobenius test (1 if null is rejected, 0 otherwise)\cr
#' tOp: output of Boot-Spectral test (1 if null is rejected, 0 otherwise)\cr
#' pvalFro: p-value for Boot-Frobenius test\cr
#' pvaOP: p-value for Boot-Spectral test

#' @export
#'
#' @examples
#' A <- genSparseGraph(4, model=list(name='ER',n=10,p=.1))
#' B <- genSparseGraph(4, model=list(name='ER',n=10,p=.5))
#' test_result = ShufflingTests(A, B, 0.05, 200)
#'
ShufflingTests <- function(A, B, sig, bs) {
  m <- min(length(A), length(B))
  C <- c(A, B)

  testStat <- Boot_Spectral_computeStat(A[1:m], B[1:m])

  # Bootstrapping via randomly permuting all labels
  bsStat <- matrix(0, nrow = bs, ncol = 2)
  for (b in 1:bs) {
    ind <- sample(length(C))
    bsStat[b, ] <- Boot_Spectral_computeStat(C[ind[1:m]], C[ind[(m+1):(2*m)]])
  }
  bsStat <- apply(bsStat, 2, sort, decreasing = TRUE)

  # p-values and acceptance/rejection
  # +0.5 for continuity correction
  pvalFro <- (sum(bsStat[, 1] >= testStat[1]) + 0.5) / bs
  pvalOp <- (sum(bsStat[, 2] >= testStat[2]) + 0.5) / bs
  tFro <- (pvalFro <= sig)
  tOp <- (pvalOp <= sig)

  return(list(tFro = tFro, tOp = tOp, pvalFro = pvalFro, pvalOp = pvalOp))
}




# Auxiliary functions


# Boot-Spectral Compute Statistics ----------------------------------------

#' Boot-Spectral Compute Statistics
#'
#' This function computes the test statistics for Boot-Spectral test, serving as an auxiliary functions for the function: ShufflingTests
#'
#' @param A1 a sparse adjacency matrix
#' @param B1 a sparse adjacency matrix
#'
#' @return \eqn{\Tau_fro} & \eqn{\Tau_spec}
#' @export
#'
#' @examples
#' A <- genSparseGraph(4, model=list(name='ER',n=10,p=.1))
#' B <- genSparseGraph(4, model=list(name='ER',n=10,p=.5))
#' testStat <- Boot_Spectral_computeStat(A, B)
#'
#'
#'
Boot_Spectral_computeStat <- function(A1, B1) {
  m1 <- floor(length(A1) / 2)
  n <- dim(A1[[1]])[1]

  SA1 <- matrix(0, n, n)
  SA2 <- matrix(0, n, n)
  SB1 <- matrix(0, n, n)
  SB2 <- matrix(0, n, n)

  for (i in 1:m1) {
    SA1 <- SA1 + A1[[i]]
    SA2 <- SA2 + A1[[i+m1]]
    SB1 <- SB1 + B1[[i]]
    SB2 <- SB2 + B1[[i+m1]]
  }

  nummat <- pracma::triu(SA1 - SB1) * pracma::triu(SA2 - SB2)
  denmat <- pracma::triu(SA1 + SB1) * pracma::triu(SA2 + SB2)
  Stat_Op = sum(nummat) / sqrt(sum(denmat))
  Stat_Fro = svd(SA1 + SA2 - SB1 - SB2)$d[1]/sqrt(norm(SA1 + SA2 + SB1 + SB2, type = "i"))
  return(c(Stat_Fro, Stat_Op))
}
