# mesoscale network testing

# required packages:
# softImpute
# irlba
# stats

# helper functions:

# logit function
Logit <- function(x){
  log((x/(1-x)))
}

# rank-truncated psuedoinverse
Tpinv <- function(M,r){
  temp <- irlba::irlba(M,r)
  dinv <- 1/temp$d
  out <- t(temp$u %*% (dinv * t(temp$v)))
  return(out)
}

# onestep projection estimator
# takes the data, dimension, test_set, masked_set (default empty)
# both sets as 'rectangles', list of (row,col)
# returns left and righthand side o/n bases
Subspace_onestep <- function(A1bar,A2bar,d,
                             test_indices,
                             masked_indices=NULL){
  n <- nrow(A1bar)
  s_ind <- matrix(0,n,n)
  s_ind[rbind(test_indices,masked_indices)] <- 1
  mrow <- which(rowSums(s_ind)>0)
  mcol <- which(colSums(s_ind)>0)
  # estimate blocks
  Chat <- A1bar[,-mcol] - A2bar[,-mcol]
  Rhat <- A1bar[-mrow,] - A2bar[-mrow,]
  Dhat <- A1bar[-mrow,-mcol] - A2bar[-mrow,-mcol]
  # estimate T
  That <- Chat %*% Tpinv(Dhat,d) %*% Rhat
  # esimate subspaces
  temp <- irlba::irlba(That,d)
  Lproj <- temp$u
  Rproj <- temp$v
  return(list(Lproj=Lproj,Rproj=Rproj))
}

# projection estimates with hard imputation
# takes the data, dimension, test_set, masked_set (default empty)
# both sets as 'rectangles', list of (row,col)
# returns left and righthand side o/n bases
Subspace_impute <- function(A1bar,A2bar,d,
                            test_indices,
                            masked_indices=NULL){
  uindices <- rbind(test_indices,masked_indices)
  # estimate blocks
  Adiff <- A1bar - A2bar
  Adiff[uindices] <- NA
  # esimate subspaces
  impdiff <- softImpute::softImpute(Adiff,rank.max=d,lambda=0,type='svd')
  # now returns full input, orthonormality enforced later
  Lproj <- impdiff$u
  Rproj <- impdiff$v
  return(list(Lproj=Lproj,Rproj=Rproj))
}

# binary_meso helper for logistic model
# takes the data, test_set, projections, variance option
# returns the pvalue and rejection decision
Binary_meso <- function(A,B,
                        test_indices,
                        test_proj,
                        var_type){
  # left/right objects
  Lproj <- test_proj$Lproj
  Rproj <- test_proj$Rproj
  # dimensions
  ss <- nrow(test_indices)
  m <- length(A)
  n <- nrow(A[[1]])
  # reset tilde if m is 1
  if(m==1){
    var_type <- 'basic'
  }
  dd <- ncol(Lproj)*ncol(Rproj)
  # orthonormal basis
  s_ind <- matrix(FALSE,n,n)
  s_ind[test_indices] <- TRUE
  UV <- (Rproj %x% Lproj)[c(s_ind),]
  # common subspace
  X <- rbind(UV,UV)
  # difference subspace
  W <- rbind(UV,-UV)
  # final design
  Xall <- cbind(X,W)
  # populate full successes
  s_vec <- rep(0,2*ss)
  for(kk in 1:m){
    s_vec <- s_vec + c(A[[kk]][test_indices],B[[kk]][test_indices])
  }
  # calculate failures
  f_vec <- m - s_vec
  # fit logistic regression
  if(var_type=='quasi'){
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::quasibinomial())
  }
  else{
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::binomial())
  }
  # get logistic regression coefficients
  gamma1 <- moda$coefficients[1:dd]
  gamma2 <- moda$coefficients[-(1:dd)]
  # estimate covariance
  # linear model predictor
  lmhat <- c(UV %*% gamma1)
  # saturated model predictor
  prop_vec <- (s_vec[1:ss] + s_vec[-(1:ss)])/(2*m)
  prop_vec_reg <- sign(prop_vec - 0.5)*pmax(abs(prop_vec - 0.5) - (1/(4*m)),0) + 0.5
  smhat <- Logit(prop_vec_reg)
  # covariances
  Ghat <- 2*(t(UV) %*% (diag(vexpit(lmhat)) %*% UV))
  Fhat <- 2*(t(UV) %*% (diag(vexpit(smhat)) %*% UV))
  # full inverse covariance
  iVhat <- Ghat %*% (solve(Fhat) %*% Ghat)
  # stat numerator
  num <- m*sum(gamma2 * (iVhat %*% gamma2))
  # overdispersion estimate as denominator
  den <- summary(moda)$dispersion
  # compile results and return
  out <- list()
  out$stat <- num/den
  out$dfs <- dd
  out$pval <- stats::pchisq(out$stat,dd,lower.tail=FALSE)
  return(out)
}


# weight_meso2 helper for Gaussian model
# takes the data, test_set, projections, variance option (basic, quasi)
# returns the pvalue and rejection decision
## avoids using rectangular properties
## assumes n x d input projection matrices
## assumes hypothesis set passed as a 2 column matrix of indices
Weight_meso <- function(A,B,
                        test_indices,
                        test_proj,
                        var_type){
  # left/right objects
  Lproj <- test_proj$Lproj
  Rproj <- test_proj$Rproj
  # dimensions
  ss <- nrow(test_indices)
  m <- length(A)
  n <- nrow(A[[1]])
  # augmented data approach
  dd <- ncol(Lproj)*ncol(Rproj)
  # orthonormal basis
  s_ind <- matrix(FALSE,n,n)
  s_ind[test_indices] <- TRUE
  UV <- svd((Rproj %x% Lproj)[c(s_ind),])$u
  # test
  X <- Y <- matrix(NA,m,dd)
  for(kk in 1:m){
    X[kk,] <- c(crossprod(UV, c(A[[kk]][test_indices])))
    Y[kk,] <- c(crossprod(UV, c(B[[kk]][test_indices]))) ### can rewrite this with sapply
  }
  # basic F test on augmented data assuming equal variance, independence
  # sse values
  sse_null <- sum((scale(rbind(X,Y),scale=F))^2)
  sse_alt <- sum((scale(X,scale=F))^2) + sum((scale(Y,scale=F))^2)
  # adjustment for tilde
  if(var_type=='quasi'){
    Xa <- (diag(ss) - tcrossprod(UV)) %*% sapply(A,function(x){x[test_indices]})
    Ya <- (diag(ss) - tcrossprod(UV)) %*% sapply(B,function(x){x[test_indices]})
    sse_alt_den <- sum(Xa^2) + sum(Ya^2)
    df_den <- 2*m*(ss - dd)
  }
  else{
    sse_alt_den <- sse_alt
    df_den <- dd*(2*(m-1))
  }
  # statistic, dfs, pvalue
  out <- list()
  out$stat <- ((sse_null - sse_alt)/dd)/(sse_alt_den/(df_den))
  out$dfs <- dd
  out$pval <- stats::pf(out$stat,df1=dd,df2=df_den,lower.tail=FALSE)
  return(out)
}

Mesoscale_test <- function(A,B,sig,
                           test_set, # a list of (row,col) # TODO: a 2-column matrix of entries
                           edge_type='weighted',# or binary
                           # other options
                           dimension,
                           # symmetric = FALSE ## TODO
                           proj_type='impute',# or onestep
                           var_type='basic', # or 'quasi'
                           masked_set=list(NULL,NULL) # as a list of (row,col) # TODO: a 2-column matrix of entries, or a list of rectangles???
){
  # parameter checking

  # does the test set specify a rectangle
  if(is.matrix(test_set)){
    rect <- FALSE
    test_indices <- test_set
  }
  else{
    rect <- TRUE
    test_indices <- as.matrix(expand.grid(test_set))
  }

  # does the masked set specify a rectangle
  if(!is.null(masked_set)){
    if(is.matrix(masked_set)){
      masked_indices <- masked_set
    }
    else{
      masked_indices <- as.matrix(expand.grid(masked_set))
      if(nrow(masked_indices)==0){
        masked_indices <- NULL
      }
    }
  }
  else{
    masked_indices <- NULL
  }

  # dimensions
  m <- length(A)
  n <- nrow(A[[1]])
  # group means
  Abar <- Reduce('+',A)/m
  Bbar <- Reduce('+',B)/m
  # transformation for binary edges
  if(edge_type=='binary'){
    delta <- 1/(3*m)
    Abar <- Logit(pmin(pmax(Abar,delta),1-delta))
    Bbar <- Logit(pmin(pmax(Bbar,delta),1-delta))
  }
  # learn left and right-hand side projections
  if(proj_type=='impute'){
    test_proj <- Subspace_impute(Abar,Bbar,dimension,test_indices,masked_indices)
  }
  else{
    test_proj <- Subspace_onestep(Abar,Bbar,dimension,test_indices,masked_indices)
  }
  # test procedures
  if(edge_type=='weighted'){
    test_out <- Weight_meso(A,B,test_indices,test_proj,var_type)
  }
  else{
    test_out <- Binary_meso(A,B,test_indices,test_proj,var_type)
  }
  # returns acceptance/rejection decision and a pvalue
  return(c(as.integer(test_out$pval <= sig),test_out$pval))
}
