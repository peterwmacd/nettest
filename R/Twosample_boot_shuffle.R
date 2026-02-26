#### helper function ####

# compute estimates of matrix norm differences
matnorm_stats <- function(A, B, op=TRUE) {
  # dimension
  m = floor(min(length(A), length(B)) / 2)

  # sum subsamples
  A1 <- Reduce('+',A[1:m])
  B1 <- Reduce('+',B[1:m])
  A2 <- Reduce('+',A[((1+m):(2*m))])
  B2 <- Reduce('+',B[((1+m):(2*m))])

  # compute statistic
  numer <- pracma::triu(A1-B1)*pracma::triu(A2-B2)
  denom <- pracma::triu(A1+B1)*pracma::triu(A2+B2)

  stat_frob = sum(numer) / sqrt(sum(denom))

  if(op){
    stat_op = irlba::irlba(A1 + A2 - B1 - B2,1)$d/sqrt(norm(A1 + A2 + B1 + B2, type = "I"))
    return(c(frob=stat_frob,op=stat_op))
  }
  else{
    return(stat_frob)
  }
}

#### main function ####

#' Two-sample Bootstrap Shuffling Tests
#'
#' Two bootstrap tests for network two-sample testing, described in \href{https://arxiv.org/abs/1811.12752}{Ghoshdastidar & von Luxburg, (2018)}.
#' Both are used to test the global null hypothesis under the undirected IER model with no self-loops.
#' The Frobenius test statistic is an estimate of the Frobenius norm of the difference in expected adjacency matrices.
#' The operator test statistic is an estimate of the operator norm of the difference in expected adjacency matrices.
#' The test is calibrated using a bootstrap which shuffles the two samples.
#' This test requires at least 2 networks in each sample.
#'
#' @param A A list of matrices, adjacency matrices for the first sample.
#' @param B A list of matrices, adjacency matrices for the second sample.
#' @param nboot number of bootstrap replicates, defaults to \code{500}.
#' @param return_null Boolean, should the full bootstrap null distribution be returned? Defaults to \code{FALSE}
#'
#' @return A list containing:
#' \item{p.value}{Named vector containing the bootstrap \eqn{p}-values for the Frobenius and operator tests (with continnuity correction)}
#' \item{statistic}{Named vector containing the test statistics for the Frobenius and operator tests}
#' \item{null_dist}{If \code{return_null=TRUE}, \code{nboot}\eqn{\times}\code{2} matrix containing the bootstrap samples for the Frobenius and operator norm tests}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=10,model=m1,m=4)$A
#' B <- Simulate_netmodel(n=10,model=m2,m=4)$A
#' test <- Twosample_boot_shuffle(A,B,nboot=40,return_null=TRUE)
Twosample_boot_shuffle <- function(A, B, nboot=500, return_null=FALSE) {
  # dimension and merge lists
  m <- min(length(A), length(B))
  C <- c(A, B)

  # observed statistics
  statistic <- matnorm_stats(A[1:m], B[1:m])

  # Bootstrapping via randomly permuting all labels
  null_dist <- matrix(0, nrow = nboot, ncol = 2)
  colnames(null_dist) <- c('frob','op')
  for (bb in seq_len(nboot)) {
    ind <- sample(length(C))
    null_dist[bb, ] <- matnorm_stats(C[ind[1:m]], C[ind[(m+1):(2*m)]])
  }
  null_dist <- apply(null_dist, 2, sort)

  # p-values
  # +0.5 for continuity correction
  p.value <- bootp_multi(statistic,null_dist)
  names(p.value) <- c('frob','op')

  # return result
  if(return_null){
    return(list(p.value=p.value,statistic=statistic,null_dist=null_dist))
  }
  else{
    return(list(p.value=p.value,statistic=statistic))
  }
}
