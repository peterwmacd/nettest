#### helper function ####

ase_stats <- function(A,B,d,return_P=FALSE) {
  # Calculate the element-wise average
  Abar = Reduce("+", A) / length(A)
  Bbar = Reduce("+", B) / length(B)

  # SVD for sample 1
  u1svd <- irlba::irlba(Abar,d)
  u1 <- u1svd$u
  s1 <- diag(u1svd$d)
  v1 <- u1svd$v

  # SVD for sample 2
  u2svd <- irlba::irlba(Bbar,d)
  u2 <- u2svd$u
  s2 <- diag(u2svd$d)
  v2 <- u2svd$v

  # compute ASEs
  X <- u1 %*% sqrt(s1)
  Y <- u2 %*% sqrt(s2)

  # ASE and EPA statistics
  statistic <- numeric(2); names(statistic) <- c('ase','epa')
  statistic[1] <- norm(X - proc_align(Y,X), type = "F")
  PA <- u1 %*% s1 %*% t(v1)
  PB <- u2 %*% s2 %*% t(v2)
  statistic[2] <- norm(PA - PB, type = "F")

  # return result
  if(return_P){
    return(list(statistic=statistic,PA=PA,PB=PB))
  }
  else{
    return(statistic)
  }
}

#### main function ####

#' Two-sample Bootstrap Tests using ASE
#'
#' Two bootstrap tests for network two-sample testing.
#' Both are used to test the global null hypothesis under the undirected IER model with no self-loops, and assuming the
#' expected adjacency matrices have rank \code{d}.
#' The ASE test statistic is based on the (Frobenius norm) difference in the aligned adjacency spectral embeddings for
#' each sample.
#' The EPA test statistic is based on the (Frobenius norm) difference in the estimated expected adjacency matrices.
#' The test is calibrated using a parametric bootstrap.
#' This test supports one network per sample.
#'
#' @param A A matrix or list of matrices, adjacency matrix(es) for the first sample.
#' @param B A matrix or list of matrices, adjacency matrix(es) for the second sample.
#' @param d rank for ASE embedding
#' @param nboot number of bootstrap replicates, defaults to \code{500}.
#' @param return_null Boolean, should the full bootstrap null distribution be returned? Defaults to \code{FALSE}
#'
#' @return A list containing:
#' \item{p.value}{Named vector containing the bootstrap \eqn{p}-values for the ASE and EPA tests (with continuity correction)}
#' \item{statistic}{Named vector containing the test statistics for the ASE and EPA tests}
#' \item{null_dist}{If \code{return_null=TRUE}, a list of 2 \code{nboot}\eqn{\times}\code{2} matrices containing the bootstrap samples for the ASE and EPA tests}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=10,model=m1)$A
#' B <- Simulate_netmodel(n=10,model=m2)$A
#' test <- Twosample_boot_ase(A,B,d=2,nboot=40,return_null=TRUE)
#'
Twosample_boot_ase <- function(A, B, d, nboot=500, return_null=FALSE) {
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  B <- checklist(B)
  mA <- length(A); mB <- length(B)

  # compute the test statistic
  ase_fit <- ase_stats(A,B,d,return_P=TRUE)
  statistic <- ase_fit$statistic

  # parametric bootstrap using estimated PA
  null_distA <- matrix(0, nrow = nboot, ncol = 2)
  colnames(null_distA) <- c('ase','epa')
  for (bb in seq_len(nboot)) {
    C <- Simulate_ier_sample(ase_fit$PA,mA+mB)
    null_distA[bb, ] <- ase_stats(C[1:mA], C[((mA+1):(mA+mB))],d)
  }
  null_distA <- apply(null_distA, 2, sort)
  # p-values using estimated PA
  p.valueA <- bootp_multi(statistic,null_distA)

  # parametric bootstrap using estimated PB
  null_distB <- matrix(0, nrow = nboot, ncol = 2)
  colnames(null_distB) <- c('ase','epa')
  for (bb in seq_len(nboot)) {
    C <- Simulate_ier_sample(ase_fit$PB,mA+mB)
    null_distB[bb, ] <- ase_stats(C[1:mA], C[((mA+1):(mA+mB))],d)
  }
  null_distB <- apply(null_distB, 2, sort)
  # p-values using estimated PA
  p.valueB <- bootp_multi(statistic,null_distB)

  # final p-values and return results
  p.value <- c(ase=max(p.valueA[1],p.valueB[1]),epa=max(p.valueA[2],p.valueB[2]))
  if(return_null){
    return(list(p.value=p.value,statistic=statistic,null_dist=list(A=null_distA,B=null_distB)))
  }
  else{
    return(list(p.value=p.value,statistic=statistic))
  }
}
