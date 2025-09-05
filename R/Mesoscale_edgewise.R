# edgewise baseline approaches to mesoscale testing
# * only implemented for rectangular hypothesis sets

# helper: convert multiplex network data (3d array) to data matrix (2d)
multnet_to_data <- function(A, # n x n x m array
                         row_indices,col_indices,
                         directed=TRUE,
                         self_loops=TRUE){
  # convert row/col to hypothesis set indices
  hyp_indices <- as.matrix(expand.grid(list(row_indices,col_indices)))
  # account for self loops in masked set and hypothesis set
  if(!self_loops){
    hyp_sl <- hyp_indices[,1]==hyp_indices[,2]
    if(sum(hyp_sl) > 0){
      hyp_indices <- hyp_indices[!hyp_sl,]
    }
  }
  # account for symmetry in masked set and hypothesis set
  if(!directed){
    # keep unique after unambigous ordering i \leq j
    hyp_indices <- unique(cbind(apply(hyp_indices,1,min),apply(hyp_indices,1,max)))
  }
  # initialize data matrix
  # dimensions
  m <- dim(A)[3]
  s <- nrow(hyp_indices)
  Z <- matrix(NA,m,s)
  # fill by row
  for(kk in 1:m){
    Z[kk,] <- A[,,kk][hyp_indices]
  }
  # return converted matrix
  return(Z)
}

# baselines

#' @export
Mesoscale_edgewise <- function(A1,A2,sig=0.05,
                          row_indices,col_indices,
                          directed=TRUE,
                          self_loops=TRUE){
  # convert lists of matrices to 3D arrays
  A1a <- list_to_array(A1)
  A2a <- list_to_array(A2)
  # convert arrays to matrix data
  Z1 <- multnet_to_data(A1a,row_indices,col_indices,directed,self_loops)
  Z2 <- multnet_to_data(A2a,row_indices,col_indices,directed,self_loops)
  # dimensions
  s <- ncol(Z1)
  m1 <- nrow(Z1)
  m2 <- nrow(Z2)
  # mean difference
  W <- colMeans(Z1) - colMeans(Z2)
  # edgewise sample variances
  V1 <- apply(Z1,2,var)
  V2 <- apply(Z2,2,var)
  # standardized statistics
  X <- W / sqrt((V1/m1) + (V2/m2))
  # combined sum statistic
  SX <- sum(X^2)
  # combined maxl statistic
  MX <- max(X^2)
  # compute test result, pval, statistic, df for each of 4 methods
  out <- list(chi=list(),Z=list(),maxl=list(),split=list())
  # Sum with chi2-approximation
  out$chi$df <- s
  out$chi$stat <- SX
  out$chi$pval <- stats::pchisq(SX,df=s,lower.tail=FALSE)
  out$chi$result <- as.integer(out$chi$pval <= sig)
  # Sum with z-approximation
  out$Z$stat <- (SX - s)/sqrt(2*s)
  out$Z$pval <- stats::pnorm(out$Z$stat,lower.tail=FALSE)
  out$Z$result <- as.integer(out$Z$pval <= sig)
  # Max'l with Xia+Lin approximation
  out$maxl$stat <- MX
  out$maxl$cutoff <- log(s) - log(log(s)) - log(pi) - 2*log(log(1/(1-sig)))
  out$maxl$pval <- 1 - exp(-exp((-1/2)*(MX - log(s) + log(log(s)) + log(pi))))
  out$maxl$result <- as.integer(MX >= out$maxl$cutoff)

  # preliminaries for sample splitting
  m11 <- floor(m1/2)
  m21 <- floor(m2/2)
  m12 <- m1 - m11
  m22 <- m2 - m21
  split1 <- sample(1:m1,floor(m1/2))
  split2 <- sample(1:m2,floor(m2/2))
  # convert arrays to matrix data
  Z11 <- multnet_to_data(A1a[,,split1],row_indices,col_indices,directed,self_loops)
  Z12 <- multnet_to_data(A1a[,,-split1],row_indices,col_indices,directed,self_loops)
  Z21 <- multnet_to_data(A2a[,,split2],row_indices,col_indices,directed,self_loops)
  Z22 <- multnet_to_data(A2a[,,-split2],row_indices,col_indices,directed,self_loops)
  # mean differences
  W1 <- colMeans(Z11) - colMeans(Z21)
  W2 <- colMeans(Z12) - colMeans(Z22)
  # combined normalizing variance
  Vsplit <- sum(((V1/m11) + (V2/m21))*((V1/m12) + (V2/m22)))
  # standardized statistic
  Ssplit <- sum(W1*W2) / sqrt(Vsplit)

  # split-based stat with normal approximation
  out$split$stat <- Ssplit
  out$split$pval <- stats::pnorm(Ssplit,lower.tail=FALSE)
  out$split$result <- as.integer(out$split$pval <= sig)

  # return all 4 test results
  return(out)
}





