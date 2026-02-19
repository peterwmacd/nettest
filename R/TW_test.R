# TW Test -----------------------------------------------------------------

# Returns:
#   acceptance/rejection decision
#   p-value

# Output:
#   A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#   B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#   sig: significance level for acceptance of null hypothesis
#   r: r communities in G

# Example: Asymp_TW(simu$A_G, simu$A_H, 0.05, 2)




#' Tracy Widom Test
#'
#'
#' The TW test is a two-sample testing procedure for large graphs.
#' The TW test is based on the Wishart distribution and uses the Tracy-Widom law to approximate the null distribution of a test statistic.\cr\cr
#' The TW test then uses the Tracy-Widom law to approximate the null distribution of T_TW under certain assumptions about the underlying distributions.
#' If the value of T_TW exceeds a certain threshold, then the null hypothesis (i.e., both populations have the same population adjacency) is rejected.
#'
#' @param A cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#' @param B cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#' @param sig significance level for acceptance of null hypothesis
#' @param r r communities in G
#'
#' @return acceptance/rejection decision & p-value for a normal distribution
#' @export
#'
#' @examples
#' A <- genSparseGraph(1,model=list(name='2SBM',n=10,p=0.8,q=0.3))
#' B <- genSparseGraph(1,model=list(name='2SBM',n=10,p=1,q=0.3))
#' test_result = Asymp_TW(A, B, 0.05, 2)
#'
Asymp_TW <- function(A, B, sig, r) {
  n <- dim(A[[1]])[1]
  m1 <- length(A)
  m2 <- length(B)

  # spectral clustering
  #idx = spectral_clus((A[[1]]+B[[1]])/2,r);

  # (temporary) estimation with graphon::est.LG
  # eventually add all options like Asymp_third_power.R
  P_1 <- graphon::est.LG(A, r)$P
  P_2 <- graphon::est.LG(B, r)$P

  # Compute sample average for the two group
  Abar = Reduce("+", A) / length(A)
  Bbar = Reduce("+", B) / length(B)
  denom = sqrt((( (P_1*(1 - P_1)) / m1) + ((P_2*(1 - P_2)) / m2)) * (n-1))
  #compute C
  C = (Abar - Bbar)/denom

  #test statistics
  spectral.norm = irlba::irlba(C, 1)$d # the largest singular value
  test.stat = (n ^ (2 / 3)) * (spectral.norm - 2)
  p.val = RMTstat::ptw(test.stat, beta=1, lower.tail = FALSE)
  test <- ifelse(p.val <= sig, 1, 0)

  return(c(test, p.val))
}

