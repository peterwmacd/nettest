# Chi-square2 Test ---------------------------------------------------------
# Note:
#   All graphs are assumed to be unweighted, undirected, and defined on a common vertex set.
#   Potential commands: pchisq(..., n * (n - 1) / 2, lower.tail = FALSE)

# Input:
#   A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#   B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#   sig: significance level for acceptance of null hypothesis

# Output:
#   acceptance/rejection decision
#   p-value for the chi2 test

# Example: Asymp_chi2(simu$A_G, simu$A_H, 0.05)

#' Chi-square2 Test
#'
#'
#' The Asymp-Chi2 test is an asymptotic test for two-sample testing of large graphs.
#' The test is used to determine whether two graphs (or graph populations) have the same population adjacency.
#' The Asymp-Chi2 test involves computing a test statistic T_chi2 (which is a sum of squared differences between the sample mean differences and their estimated variances) \cr\cr
#' The test involves computing rejection decision and the p-value, which is the probability that the null hypothesis (i.e., both graphs have the same population adjacency) is true.
#'
#' @param A cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#' @param B cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#' @param sig significance level for acceptance of null hypothesis
#' @param cov_method covariance estimation method, one of 'diag', 'shrink', or 'full'. Full covariance
#' will be singular unless there are more total samples than edges
#'
#' @return acceptance/rejection decision & p-value for a normal distribution
#' @export
#'
#' @examples
#' A <- genSparseGraph(4,model=list(name='2SBM',n=10,p=0.5,q=0.3))
#' B <- genSparseGraph(4,model=list(name='2SBM',n=10,p=0.5,q=0.3))
#' test_result = Asymp_chi2(A, B, 0.05)
#'
Asymp_chi2 <- function(A, B, sig, cov_method='diag') {

  m.a = length(A)
  m.b = length(B)
  n <- dim(A[[1]])[1] # nodes

  vec.A = array(0, dim = c(m.a, n * n))
  vec.B = array(0, dim = c(m.b, n * n))

  for (i in 1:m.a){
    curr = A[[i]]
    curr[lower.tri(curr, diag = TRUE)] = 0
    vec.A[i, ] = t(as.vector(curr))
  }
  for (i in 1:m.b){
    curr = B[[i]]
    curr[lower.tri(curr, diag = TRUE)] = 0
    vec.B[i, ] = t(as.vector(curr))
  }

  idx_non_zero = (colSums(vec.A + vec.B) != 0)
  vec.A = vec.A[, idx_non_zero]
  vec.B = vec.B[, idx_non_zero]

  mean.diff = colMeans(vec.A) - colMeans(vec.B)
  #numer = mean.diff ^ 2

  denom = apply(vec.A, MARGIN = 2, stats::var) / m.a +
   apply(vec.B, MARGIN = 2, stats::var) / m.b
  # remove zero variances for covariance estimation
  ind_nonzero_dem = (denom != 0)
  mean.diff <- mean.diff[ind_nonzero_dem]
  vec.A <- vec.A[,ind_nonzero_dem]
  vec.B <- vec.B[,ind_nonzero_dem]

  # cov_method one of diag, shrink, full
  # previous method is diag (scale by marginal variances assuming independent edges)
  if(cov_method=='diag'){
    Q <- diag(1/denom[ind_nonzero_dem])
  }
  else if(cov_method=='shrink'){
    # shrinkage estimation (Shafer + Strimmer 2005), as in Ginestet et al
    suppressWarnings(sigma.A <- corpcor::cov.shrink(vec.A,verbose=FALSE) / m.a)
    suppressWarnings(sigma.B <- corpcor::cov.shrink(vec.B,verbose=FALSE) / m.b)
    Q <- solve(sigma.A + sigma.B)
  }
  else if(cov_method=='full'){
    if(ncol(vec.A) > (m.a + m.b)){
      print("insufficient sample size for full covariance")
    }
    sigma.A <- stats::cov(vec.A) / m.a
    sigma.B <- stats::cov(vec.B) / m.b
    Q <- solve(sigma.A + sigma.B)
  }
  else{print("unexpected input for cov_method")}

  #denom = 0
  #ind_nonzero_dem = (denom != 0)
  #numer = numer[ind_nonzero_dem]
  #denom = denom[ind_nonzero_dem]

  #test.stat = sum(numer / denom)
  test.stat <- sum( mean.diff * (Q %*% mean.diff))
  p.val = stats::pchisq(test.stat, df = n * (n - 1) / 2, lower.tail = FALSE)
  test <- ifelse(p.val <= sig, 1, 0)
  return(c(test, p.val))
}
