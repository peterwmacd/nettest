#' Two-sample Asymptotic Normal Test
#'
#' An asymptotic test for network two-sample testing.
#' It is used to test the global null hypothesis under the undirected IER model with no self-loops.
#' The test is based on an estimate of the Frobenius norm of the difference in expected adjacency matrices,
#' estimated by sample splitting.
#' The test statistic is asymptotically normal as \eqn{n \rightarrow \infty}.
#' This test requires at least 2 networks in each sample.
#'
#'
#' @param A A list of matrices, adjacency matrices for the first sample.
#' @param B A list of matrices, adjacency matrices for the second sample.
#'
#' @return A list containing:
#' \item{p.value}{The (two-sided) \eqn{p}-value for the test}
#' \item{statistic}{The asymptotically normal test statistic}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=10,model=m1,m=4)$A
#' B <- Simulate_netmodel(n=10,model=m2,m=4)$A
#' test <- Twosample_anorm(A,B)
Twosample_anorm <- function(A, B) {
  # compute statistic (using matnorm_stats from Twosample_boot_shuffle.R)
  statistic <- matnorm_stats(A,B,op=FALSE)
  p.value <- 2*stats::pnorm(statistic, lower.tail = FALSE)

  # return result
  return(list(p.value=p.value,
              statistic=statistic))
}
