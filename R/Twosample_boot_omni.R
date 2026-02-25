#### helper functions ####

# construct omni matrix
make_omni_matrix <- function(C){
  # dimensions
  n <- nrow(C[[1]])
  m <- length(C)
  # initialize omni matrix
  M <- matrix(0,n*m,n*m)
  # fill matrix
  for(ii in 1:m){
    rstart <- n*(ii-1)+1; rend <- n*ii
    for(jj in 1:m){
      cstart <- n*(jj-1)+1; cend <- n*jj
      # populate
      M[rstart:rend,cstart:cend] <- 0.5*(C[[ii]] + C[[jj]])
    }
  }
  return(M)
}

# compute statistic from omni embedding
omni_stat <- function(C,mA,mB,d) {
  # dimensions and concatenate samples
  n <- nrow(C[[1]])
  # Matrix M
  M <- make_omni_matrix(C)

  # Obtain omni embedding
  Msvd <- irlba::irlba(M,d)
  XM <- Msvd$u %*% diag(sqrt(Msvd$d))

  # get group-averages of omni embeddings
  rbreak <- n*mA
  XA <- avg_block(XM[1:rbreak,],mA)
  XB <- avg_block(XM[(1+rbreak):nrow(XM),],mB)

  # ASE statistic
  statistic <- norm(XA - XB, type = "F")
  return(statistic)
}

#### main function ####

#' Two-sample Bootstrap Tests using omnibus embedding
#'
#' Bootstrap test for network two-sample testing.
#' Used to test the global null hypothesis under the undirected IER model with no self-loops, and assuming the
#' expected adjacency matrices have rank \code{d}.
#' The ASE test statistic is based on the (Frobenius norm) difference in the omnibus (OMNI) embeddings for
#' each sample.
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
#' \item{p.value}{The bootstrap \eqn{p}-value for the OMNI test (with continuity correction)}
#' \item{statistic}{Named vector containing the test statistics for the ASE and EPA tests}
#' \item{null_dist}{If \code{return_null=TRUE}, a list of 2 \code{nboot}-vectors containing the bootstrap samples for the OMNI tests}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=10,model=m1)$A
#' B <- Simulate_netmodel(n=10,model=m2)$A
#' test <- Twosample_boot_omni(A,B,d=2,nboot=40,return_null=TRUE)
#'
Twosample_boot_omni <- function(A, B, d, nboot=500, return_null=FALSE) {
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  B <- checklist(B)
  mA <- length(A); mB <- length(B)
  C <- c(A,B)

  # compute test statistic
  statistic <- omni_stat(C,mA,mB,d)

  # ASE fitting for bootstrap calibration
  ase_fit <- ase_stats(A,B,d,return_P=TRUE)

  # parametric bootstrap using estimated PA
  null_distA <- numeric(nboot)
  for (bb in seq_len(nboot)) {
    Cboot <- Simulate_ier_sample(ase_fit$PA,mA+mB)
    null_distA[bb] <- omni_stat(Cboot,mA,mB,d)
  }
  null_distA <- sort(null_distA)
  # p-values using estimated PA
  p.valueA <- bootp(statistic,null_distA)

  # parametric bootstrap using estimated PB
  null_distB <- numeric(nboot)
  for (bb in seq_len(nboot)) {
    Cboot <- Simulate_ier_sample(ase_fit$PB,mA+mB)
    null_distB[bb] <- omni_stat(Cboot,mA,mB,d)
  }
  null_distB <- sort(null_distB)
  # p-values using estimated PA
  p.valueB <- bootp(statistic,null_distB)

  # final p-values and return results
  p.value <- max(p.valueA,p.valueB)
  if(return_null){
    return(list(p.value=p.value,statistic=statistic,null_dist=list(A=null_distA,B=null_distB)))
  }
  else{
    return(list(p.value=p.value,statistic=statistic))
  }
}




