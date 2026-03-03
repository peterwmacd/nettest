#' Two-sample 3rd Power Test
#'
#' An asymptotic test for network two-sample testing, described in \href{https://www.jstor.org/stable/27280115}{Chen et al., (2024)}.
#' It is used to test the global null hypothesis under the undirected IER model with no self-loops.
#' The test is based on the trace of the 3rd power of a suitably scaled difference of mean adjacency matrices, with a random adjustment to the diagonal entries.
#' The test statistic is asymptotically normal as \eqn{n \rightarrow \infty}, as long as \eqn{m} does not grow too fast.
#'
#' @param A A matrix or list of matrices, adjacency matrix(es) for the first sample.
#' @param B A matrix or list of matrices, adjacency matrix(es) for the second sample.
#' @param P_method P-matrix estimation method, one of \code{'NBS'}, \code{'SBA'}, \code{'USVT'} or \code{spectral_clust}.
#' All but \code{spectral_clust} are implemented using the \code{graphon} package.
#' Defaults to \code{'SBA'}. Note that \code{'NBS'} and \code{'USVT'} can perform erratically when \eqn{n} is very small.
#' @param delta precision parameter larger than zero, only used if \code{P_method='SBA'}. Defaults to \code{0.5}.
#' @param d number of clusters, only used if \code{P_method='spectral_clust'}. Defaults to \code{2}.
#' @param ndiag number of random replications of the diagonal adjustment. Defaults to \code{100}.
#'
#' @return A list containing:
#' \item{p.value}{The (two-sided) \eqn{p}-value for the test}
#' \item{statistic}{The asymptotically normal test statistic, averaged over random replications of the diagonal adjustment.}
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=40,model=m1,m=4)$A
#' B <- Simulate_netmodel(n=40,model=m2,m=4)$A
#' test <- Twosample_thirdpower(A,B,P_method='SBA')
Twosample_thirdpower <- function(A, B, P_method='SBA',delta=0.5,d=2,ndiag=100){
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  B <- checklist(B)
  n <- nrow(A[[1]])
  mA <- length(A); mB <- length(B)

  # compute group sample averages
  Abar = Reduce("+", A) / mA
  Bbar = Reduce("+", B) / mB

  # P_estimation for groups 1 and 2
  switch(P_method,
         SBA={
           # NOTE: for m > 1, SBA returns P \in [0,m]
           PA = graphon::est.SBA(A, delta)$P / mA
           PB = graphon::est.SBA(B, delta)$P / mB
         },
         NBD={
           # NOTE: for m > 1, nbdsmooth returns P \in [0,1]
           PA = graphon::est.nbdsmooth(A)$P
           PB = graphon::est.nbdsmooth(B)$P
         },
         USVT={
           # NOTE: for m > 1, SBA returns P \in [0,1]
           PA = graphon::est.USVT(A, eta=0.001)$P
           PB = graphon::est.USVT(B, eta=0.001)$P
         },
         spectral_clust={
           # cluster on pooled mean adjacency matrix
           C <- spectral_clust((Abar+Bbar)/2,d)
           Z <- C_to_Z(C,d)
           # store estimates
           PA=block_avg(Abar,Z)
           PB=block_avg(Bbar,Z)
         },
         {
           stop('unexpected input for P_method')
         })

  # compute group sample averages
  Abar = Reduce("+", A) / mA
  Bbar = Reduce("+", B) / mB

  # compute normalized adjacency difference
  numer = Abar - Bbar
  denom = sqrt((( (PA*(1 - PA)) / mA) + ((PB*(1 - PB)) / mB)) * n)
  # correction to avoid zero variance for sparse entries
  denom[numer==0] <- 1
  Z = numer / denom

  # populate diagonal entries of Z for random replications
  statvec <- numeric(ndiag)
  for(ii in seq_len(ndiag)){
    Z_temp <- Z
    diag(Z_temp) <- (1 - 2*stats::rbinom(n,1,0.5))/sqrt(n)
    Z_cubed <- (Z_temp %*% Z_temp) %*% Z_temp
    statvec[ii] <- sum(diag(Z_cubed)) / sqrt(15)
  }

  # summarize statistic
  statistic <- mean(statvec)
  p.value <- 2*stats::pnorm(abs(statistic),lower.tail=FALSE)

  # return result
  return(list(p.value=p.value,
              statistic=statistic))
}

