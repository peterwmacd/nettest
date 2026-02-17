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

# 4 edgewise approaches (see response letter)
# (1 chi) aggregate (sum) edgewise chi-squared stats
# (2 Z) aggregate (sum) edgewise chi-squared stats with normal approx as rc -> \infty
# (3 maxl) aggregate (max) edgewise chi-squared stats with extreme value thresh
# (4 split) aggregate (sum) products of mean differences, sample split to avoid 4th moment
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

# 1D projection based on A1-A2, plus recalibration for double-dipping bias
# NOTE only implemented for directed networks with self loops
#' @export
Mesoscale_1d <- function(A1,A2,
                         sig=0.05,nboot=1e3,
                         sigma=NULL,
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
  # mean difference (normalized to a direction)
  W <- colMeans(Z1) - colMeans(Z2)
  W <- W / sqrt(sum(W^2))
  # pooled statistic (squared discrepancy, equivalent to 2-sided normal test)
  D <- (mean((Z1 - Z2) %*% W))^2

  # estimated sample variance from entrywise variances
  if(is.null(sigma)){
    sigma_use <- mean(c(apply(Z1,2,var),apply(Z2,2,var)))
  }
  else{
    sigma_use <- sigma
  }

  # theoretical (selection biased cutoff)
  # null variance of mean(Z1*W) is sigma^2/m1; variance of mean(Z2*W) is sigma^2/m2
  # null variance of D is (sigma^2)*(1/m1 + 1/m2)
  # cutoff is based on a one-sided test, upper 1-sig quantile of normal(0,var) = sd*qnorm(1-sig)
  sigma_theory <- sigma_use*((1/m1) + (1/m2))

  # bootstrap calibration
  DB <- rep(NA,nboot)
  for(bb in 1:nboot){
    Z1B <- matrix(rnorm(s*m1,sd=sqrt(sigma_use)),ncol=s)
    Z2B <- matrix(rnorm(s*m2,sd=sqrt(sigma_use)),ncol=s)
    WB <- colMeans(Z1B) - colMeans(Z2B)
    WB <- WB / sqrt(sum(WB^2))
    # pooled statistic
    DB[bb] <- (mean((Z1B - Z2B) %*% WB))^2
  }

  # theoretical and bootstrap calibrated cutoff and pval
  out <- list(theory=list(),boot=list())

  # theoretical results
  out$theory$stat <- D
  out$theory$cutoff <- sigma_theory*qchisq(sig,df=1,lower.tail=FALSE)
  out$theory$pval <- pchisq(D/sigma_theory,df=1,lower.tail=FALSE)
  out$theory$result <- as.integer(D >= out$theory$cutoff)

  # bootstrap results
  out$boot$stat <- D
  out$boot$cutoff <- as.numeric(quantile(DB,1-sig))
  out$boot$pval <- mean(DB > D)
  out$boot$result <- as.integer(D >= out$boot$cutoff)

  # return test results
  return(out)
}

# two directions based on data

# oracle projection using Theta1, Theta2
orc_dir <- function(Theta1,Theta2,row_indices,col_indices){
  v <- c(Theta1[row_indices,col_indices] - Theta2[row_indices,col_indices])
  v <- v/sqrt(sum(v^2))
  return(v)
}

# imputation-based direction using mean(A1),mean(A2)
# requires an additional dimension (for imputation)
impute_dir <- function(A1,A2,row_indices,col_indices,d){
  # convert lists of matrices to 3D arrays
  A1a <- list_to_array(A1)
  A2a <- list_to_array(A2)
  # compute means of arrays
  A1bar <- apply(A1a,c(1,2),mean)
  A2bar <- apply(A2a,c(1,2),mean)
  Adiff <- A1bar - A2bar
  # fill in NAs for missing entries
  Adiff[row_indices,col_indices] <- NA
  # estimate subspaces
  impdiff <- softImpute::softImpute(Adiff,rank.max=d,lambda=0,type='svd')
  # reconstruct matrix
  Thetadiff_hat <- impdiff$u %*% tcrossprod(diag(impdiff$d),impdiff$v)
  # recover direction SOLN 1:
  v <- c(Thetadiff_hat[row_indices,col_indices])
  v <- v/sqrt(sum(v^2))
  return(v)
}

# 1D with projection provided ahead; note this approach is only adapted to
# directed networks with self loops
#' @export
Mesoscale_1d_fix <- function(A1,A2,
                         sig=0.05,proj_dir,
                         sigma=NULL,
                         row_indices,col_indices){
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
  # pooled statistic
  D <- (mean((Z1 - Z2) %*% proj_dir))^2

  # estimated sample variance from entrywise variance, or use true variance
  # if provided
  if(is.null(sigma)){
    sigma_use <- mean(c(apply(Z1,2,var),apply(Z2,2,var)))
  }
  else{
    sigma_use <- sigma
  }

  # theoretical (selection biased cutoff)
  # null variance of mean(Z1*W) is sigma^2/m1; variance of mean(Z2*W) is sigma^2/m2
  # null variance of D is (sigma^2)*(1/m1 + 1/m2)
  # cutoff is based on a one-sided test, upper 1-sig quantile of normal(0,var) = sd*qnorm(1-sig)
  sigma_theory <- sigma_use*((1/m1) + (1/m2))

  # theoretically calibrated cutoff and pval
  out <- list()

  # theoretical results
  out$stat <- D
  out$cutoff <- sigma_theory*qchisq(sig,df=1,lower.tail=FALSE) # two-sided
  out$pval <- pchisq(D/sigma_theory,df=1,lower.tail=FALSE)
  out$result <- as.integer(D >= out$cutoff)

  # return test results
  return(out)
}
