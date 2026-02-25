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
#' An asymptotic test for network two-sample testing.
#' It is used to test the global null hypothesis under the undirected IER model with no self-loops.
#' The test is based on the largest singular value of a suitably scaled difference of mean adjacency matrices.
#' The test statistic is asymptotically Tracy-Widom distributed with parameter \eqn{\beta=1} as \eqn{n \rightarrow \infty}, as long as \eqn{m} does not grow too fast.
#'
#' @param A A matrix or list of matrices, adjacency matrix(es) for the first sample.
#' @param B A matrix or list of matrices, adjacency matrix(es) for the second sample.
#' @param sig significance level for acceptance of null hypothesis
#' @param r r communities in G
#'
#' @return acceptance/rejection decision & p-value for a normal distribution
#'
#' @export
#'
#' @examples
#' m1 <- list(name='SBM',B=matrix(0.3,2,2),Pi=c(0.5,0.5))
#' m2 <- m1; m2$B <- matrix(c(1,0.3,0.3,1),2,2)
#' A <- Simulate_netmodel(n=40,model=m1,m=4)$A
#' B <- Simulate_netmodel(n=40,model=m2,m=4)$A
#' test <- Twosample_tw(A,B,P_method='SBA')
TW_test <- function(A, B, P_method='SBA',delta=0.2){
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  B <- checklist(B)
  n <- nrow(A[[1]])
  mA <- length(A); mB <- length(B)

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
         {
           stop('unexpected input for P_method')
         })

  # compute group sample averages
  Abar = Reduce("+", A) / mA
  Bbar = Reduce("+", B) / mB

  # compute normalized adjacency difference
  numer = Abar - Bbar
  denom = sqrt((( (PA*(1 - PA)) / mA) + ((PB*(1 - PB)) / mB)) * (n-1))
  # correction to avoid zero variance for sparse entries
  denom[numer==0] <- 1
  C = numer / denom

  # test statistic and result
  opnorm = irlba::irlba(C, 1)$d
  statistic = (n ^ (2 / 3)) * (opnorm - 2)
  p.value = RMTstat::ptw(test.stat, beta=1, lower.tail = FALSE)

  # return result
  return(list(p.value=p.value,
              statistic=statistic))
}

