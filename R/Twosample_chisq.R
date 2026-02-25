#' Two-sample Chi-squared Test
#'
#' An asymptotic test for network two-sample testing.
#' It is used to test the global null hypothesis under the undirected IER model with no self-loops.
#' The test is based on the sum of squared difference in edgewise means.
#' The test statistic is asymptotically chi-squared distributed for fixed \eqn{n}, as \eqn{m \rightarrow \infty} (where \eqn{m} is the number of networks in each sample).
#' This test requires at least 2 networks in each sample; it will tend to perform poorly for small \eqn{m} or sparse graphs.
#'
#' @param A A list of matrices, adjacency matrices for the first sample.
#' @param B A list of matrices, adjacency matrices for the second sample.
#' @param cov_method covariance estimation method, one of \code{'diag'}, \code{'shrink'}, or \code{'full'}. Full covariance
#' will be singular unless there are more total samples than edges. Defaults to \code{'diag'}.
#'
#' @return A list containing:
#' \item{p.value}{The \eqn{p}-value for the test}
#' \item{statistic}{The asymptotically chi-square test statistic}
#' \item{df}{The degrees of freedom of the chi-square null distribution}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=10,model=m1,m=4)$A
#' B <- Simulate_netmodel(n=10,model=m2,m=4)$A
#' test <- Twosample_chisq(A,B,cov_method='shrink')
Twosample_chisq <- function(A, B, cov_method='diag') {
  # dimension
  mA <- length(A)
  mB <- length(B)
  n <- nrow(A[[1]])

  # vectorized adjacency matrices
  vecA <- list_to_vecdata(A)
  vecB <- list_to_vecdata(B)

  # remove zero variance columns for covariance estimation
  denom <- apply(vecA, 2, stats::var) / mA + apply(vecB, 2, stats::var) / mB
  ind_nonzero_dem <- (denom != 0)
  vecA <- vecA[,ind_nonzero_dem]
  vecB <- vecB[,ind_nonzero_dem]
  mean.diff <- colMeans(vecA) - colMeans(vecB)

  # cov_method one of diag, shrink, full
  switch(cov_method,
         diag = {
           Q <- diag(1/denom[ind_nonzero_dem])
         },
         shrink = {
           # shrinkage estimation (Shafer + Strimmer 2005), as in Ginestet et al
           suppressWarnings(SigmaA <- corpcor::cov.shrink(vecA,verbose=FALSE) / mA)
           suppressWarnings(SigmaB <- corpcor::cov.shrink(vecB,verbose=FALSE) / mB)
           Q <- solve(SigmaA + SigmaB)
         },
         full = {
           # check dimension for singular covariance
           if(ncol(vecA) > mA+mB){
             stop('Insufficient samples for full covariance estimation')
           }
           SigmaA <- stats::cov(vecA) / mA
           SigmaB <- stats::cov(vecB) / mB
           Q <- solve(SigmaA + SigmaB)
         },
         {
           stop('Unexpected input for cov_method')
         })

  # compute test statistic and p-value
  statistic <- sum(mean.diff * (Q %*% mean.diff))
  p.value <- stats::pchisq(statistic, df = choose(n,2), lower.tail = FALSE)
  df <- choose(n,2)
  return(list(p.value=p.value,statistic=statistic,df=df))
}
