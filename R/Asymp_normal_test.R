# Normality Test ----------------------------------------------------------

# Note:
#   All graphs are assumed to be unweighted, undirected, and defined on a common vertex set.
#   Sample size >= 2
#   Potential commands: p-value = 2 * pnorm(..., lower.tail = FALSE)

# Input:
#   A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#   B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#   sig: significance level for acceptance of null hypothesis

# Output:
#   acceptance/rejection decision & p-value for a normal distribution

# Example: Asymp_normal(simu$A_G, simu$A_H, 0.05)





#' Asymptotic Normal Test
#'
#'
#' The Asymp-Normal test is an asymptotic test for two-sample testing of large graphs.
#' It is used to determine whether two graphs (or graph populations) have the same population adjacency.\cr\cr
#' The test involves computing rejection decision and the p-value, which is the probability that the null hypothesis (i.e., both graphs have the same population adjacency) is true.
#'
#' @param A cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#' @param B cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#' @param sig significance level for acceptance of null hypothesis
#'
#' @return acceptance/rejection decision & p-value for a normal distribution
#' @export
#'
#' @examples
#' A <- genSparseGraph(4,model=list(name='2SBM',n=4,p=0.3,q=0.3))
#' B <- genSparseGraph(4,model=list(name='2SBM',n=4,p=1,q=0.3))
#' test_result = Asymp_normal(A,B, 0.05)
Asymp_normal <- function(A, B, sig) {
  m_1 = floor(min(length(A), length(B)) / 2) # m/2 (int)
  n <- dim(A[[1]])[1] # nodes

  A1 <- Matrix::Diagonal(n, 0)
  B1 <- Matrix::Diagonal(n, 0)
  A2 <- Matrix::Diagonal(n, 0)
  B2 <- Matrix::Diagonal(n, 0)

  for (i in 1:m_1) {
    #sum first part  (k <= m/2)
    A1 <- A1 + A[[i]]
    B1 <- B1 + B[[i]]

    #sum second part  (k > m/2)
    A2 <- A2 + A[[i + m_1]]
    B2 <- B2 + B[[i + m_1]]
  }

  numer = (A1 - B1)

  #i<j upper
  num.P1 = A1 - B1
  num.P2 = A2 - B2
  num.P1[lower.tri(num.P1, diag = TRUE)] = 0
  num.P2[lower.tri(num.P2, diag = TRUE)] = 0
  numer = num.P1 * num.P2

  dem.P1 = A1 + B1
  dem.P2 = A2 + B2
  dem.P1[lower.tri(dem.P1, diag = TRUE)] = 0
  dem.P2[lower.tri(dem.P2, diag = TRUE)] = 0
  denom = dem.P1 * dem.P2
  # print(sum(numer))
  # print(sum(denom))
  # print("\n")
  test.stat <- sum(numer)/sqrt(sum(denom))

  p.val <- 2 * stats::pnorm(test.stat, lower.tail = FALSE)
  test <- ifelse(p.val <= sig, 1, 0) # 1: reject, 0:not reject
  return(c(test, p.val))
}
