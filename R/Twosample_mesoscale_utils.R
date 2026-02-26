# helper functions for two-sample mesoscale testing
# - hypothesis set processing
# - one-step subspace estimation
# - imputation-based subspace estimation (using softImpute)
# - Binary edge mesoscale testing (internal)
# - Weighted edge mesoscale testing (internal)

# helper to convert different formats of hypothesis and masked sets
# of node pairs to unified 2-column matrix of (i,j)'s
# accounts for directedness and self loops in original network data
hyp_set_to_indices <- function(hyp_set,masked_set,directed,self_loops){
  # does the test set specify a rectangle, if so expand to general indices
  if(is.matrix(hyp_set)){
    hyp_indices <- hyp_set
  }
  else{
    # else list of rectangles
    if(is.list(hyp_set[[1]])){
      hyp_indices <- Reduce(rbind,lapply(hyp_set,function(x){as.matrix(expand.grid(x))}))
    }
    else{
      # else one rectangle
      hyp_indices <- as.matrix(expand.grid(hyp_set))
    }
  }
  # does the masked set specify a rectangle, if so expand to general indices
  if(!is.null(masked_set)){
    if(is.matrix(masked_set)){
      masked_indices <- masked_set
    }
    else{
      # else list of rectangles
      if(is.list(masked_set[[1]])){
        masked_indices <- Reduce(rbind,lapply(masked_set,function(x){as.matrix(expand.grid(x))}))
      }
      else{
        # else one rectangle
        masked_indices <- as.matrix(expand.grid(masked_set))
        if(nrow(masked_indices)==0){
          masked_indices <- NULL
        }
      }
    }
  }
  else{
    masked_indices <- NULL
  }
  # account for self loops in masked set and hypothesis set (add diagonal to masked set)
  if(!self_loops){
    hyp_sl <- hyp_indices[,1]==hyp_indices[,2]
    if(sum(hyp_sl) > 0){
      masked_indices <- rbind(masked_indices,hyp_indices[hyp_sl,])
      hyp_indices <- hyp_indices[!hyp_sl,]
    }
  }
  # account for symmetry in masked set and hypothesis set
  if(!directed){
    # augment hypothesis indices and masked indices to include diagonal mirrors
    hyp_indices <- unique(rbind(hyp_indices,hyp_indices[,c(2,1)]))
    if(!is.null(masked_indices)){
      masked_indices <- unique(rbind(masked_indices,masked_indices[,c(2,1)]))
    }
  }
  return(list(hyp_indices=hyp_indices,masked_indices=masked_indices))
}

# onestep projection estimator
# takes the data, dimension, hyp_indices, masked_indices (default empty)
# returns left and righthand side o/n bases
Subspace_onestep <- function(A1bar,A2bar,d,
                             hyp_indices,
                             masked_indices,
                             self_loops,
                             directed,
                             centered){
  # if no self loops, check for diagonal NAs and replace with 0's
  if(!self_loops & any(is.na(c(diag(A1bar),diag(A2bar))))){
    diag(A1bar) <- 0
    diag(A2bar) <- 0
  }
  # difference matrix
  Adiff <- A1bar - A2bar
  # dimension
  n <- nrow(A1bar)
  # compute masked rows/columns
  s_ind <- matrix(0,n,n)
  s_ind[rbind(hyp_indices,masked_indices)] <- 1
  mrow <- which(rowSums(s_ind)>0)
  mcol <- which(colSums(s_ind)>0)
  # estimate blocks
  Chat <- Adiff[,-mcol]
  Rhat <- Adiff[-mrow,]
  Dhat <- Adiff[-mrow,-mcol]
  # overall mean
  if(centered){
    mu <- mean(c(Chat,Rhat,Dhat))
  }
  else{
    mu <- 0
  }
  # estimate T
  That <- (Chat - mu) %*% Tpinv((Dhat - mu),d) %*% (Rhat - mu)
  # estimate subspaces
  temp <- irlba::irlba(That,d)
  # store left and right-hand side projections
  if(!directed){
    Lproj <- Rproj <- temp$u
  }
  else{
    Lproj <- temp$u
    Rproj <- temp$v
  }
  return(list(Lproj=Lproj,Rproj=Rproj))
}

# projection estimates with hard imputation
# takes the data, dimension, hyp_indices, masked_indices (default empty)
# returns left and righthand side o/n bases
Subspace_impute <- function(A1bar,A2bar,d,
                            hyp_indices,
                            masked_indices,
                            self_loops,
                            directed,
                            centered){
  uindices <- rbind(hyp_indices,masked_indices)
  # estimate blocks
  Adiff <- A1bar - A2bar
  # mask self loops, hypothesis set and masked indices
  if(!self_loops){
    diag(Adiff) <- NA
  }
  Adiff[uindices] <- NA
  # centering
  if(centered){
    Adiff <- Adiff - mean(Adiff,na.rm=TRUE)
  }
  # estimate subspaces
  impdiff <- softImpute::softImpute(Adiff,rank.max=d,lambda=0,type='svd')
  # store left and right-hand side projections
  if(!directed){
    Lproj <- Rproj <- impdiff$u
  }
  else{
    Lproj <- impdiff$u
    Rproj <- impdiff$v
  }
  return(list(Lproj=matrix(Lproj,ncol=d),Rproj=matrix(Rproj,ncol=d)))
}

# binary_meso helper for logistic model
# takes the data, hyp_set, projections, variance option
# returns the pvalue and rejection decision
Binary_meso <- function(A,B,
                        hyp_indices,
                        hyp_proj,
                        var_type,
                        directed,
                        centered){
  # left/right objects
  Lproj <- hyp_proj$Lproj
  Rproj <- hyp_proj$Rproj
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # reset tilde if m is 1
  if(min(m1,m2)==1){
    var_type <- 'basic'
  }
  # hypothesis set entries as a matrix
  s_ind <- matrix(FALSE,n,n)
  s_ind[hyp_indices] <- TRUE
  # 'basic' design
  if(!directed){
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind & upper.tri(s_ind,diag=TRUE)
    # linear mapping and projection basis
    Gd <- Sym_span(s_ind,s_ind_data)
    UV <- pracma::orth(crossprod(Gd,(Rproj %x% Lproj)[c(s_ind),]))
    # need orth to remove possible rank deficiency
  }
  else{
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind
    # projection basis
    UV <- (Rproj %x% Lproj)[c(s_ind),,drop=FALSE]
  }
  # centering
  if(centered){
    UV <- cbind(rep(1/sqrt(nrow(UV)),nrow(UV)),center_orth(UV))
  }
  # dimensions
  ss <- nrow(UV)
  dd <- ncol(UV)
  # common subspace
  X <- rbind(UV,UV)
  # difference subspace
  W <- rbind(UV,-UV)
  # final design
  Xall <- cbind(X,W)

  # populate successes
  s_vecA <- rowSums(sapply(1:m1,function(kk){c(A[[kk]][s_ind_data],rep(0,ss))}))
  s_vecB <- rowSums(sapply(1:m2,function(kk){c(rep(0,ss),B[[kk]][s_ind_data])}))
  s_vec <- s_vecA + s_vecB
  # populate failures
  f_vec <- c(rep(m1,ss),rep(m2,ss)) - s_vec

  # fit logistic regression
  if(var_type=='quasi'){
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::quasibinomial())
  }
  else{
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::binomial())
  }
  # get logistic regression coefficients
  gamma1 <- moda$coefficients[1:dd]
  if(centered){
    gamma2 <- moda$coefficients[-(1:(dd+1))]
  }
  else{
    gamma2 <- moda$coefficients[-(1:dd)]
  }
  # estimate covariance
  # linear model predictor
  lmhat <- c(UV %*% gamma1)
  # saturated model predictor
  prop_vec_reg <- (s_vec[1:ss] + s_vec[-(1:ss)] + 1)/(m1 + m2 + 2)
  smhat <- logit(prop_vec_reg)
  # covariances
  if(centered){
    Ghat <- 2*(t(UV[,-1,drop=FALSE]) %*% (diag(vexpit(lmhat)) %*% UV[,-1,drop=FALSE]))
    Fhat <- 2*(t(UV[,-1,drop=FALSE]) %*% (diag(vexpit(smhat)) %*% UV[,-1,drop=FALSE]))
  }
  else{
    Ghat <- 2*(t(UV) %*% (diag(vexpit(lmhat)) %*% UV))
    Fhat <- 2*(t(UV) %*% (diag(vexpit(smhat)) %*% UV))
  }
  # full inverse covariance
  iVhat <- Ghat %*% (solve(Fhat) %*% Ghat)
  # stat numerator
  num <- ((m1+m2)/2)*sum(gamma2 * (iVhat %*% gamma2))
  # overdispersion estimate as denominator
  den <- summary(moda)$dispersion
  # compile results and return
  out <- list()
  out$stat <- num/den
  if(centered){
    out$dfs <- dd-1
    out$pval <- stats::pchisq(out$stat,dd-1,lower.tail=FALSE)
  }
  else{
    out$dfs <- dd
    out$pval <- stats::pchisq(out$stat,dd,lower.tail=FALSE)
  }
  return(out)
}

# weight_meso2 helper for Gaussian model
# takes the data, hyp_set, projections, variance option (basic, quasi)
# returns the pvalue and rejection decision
## avoids using rectangular properties
## assumes n x d input projection matrices
## assumes hypothesis set passed as a 2 column matrix of indices
Weight_meso <- function(A,B,
                        hyp_indices,
                        hyp_proj,
                        var_type,
                        directed,
                        centered){
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # left/right objects
  Lproj <- hyp_proj$Lproj
  Rproj <- hyp_proj$Rproj
  # hypothesis set indices as a matrix
  s_ind <- matrix(FALSE,n,n)
  s_ind[hyp_indices] <- TRUE
  # projection matrix, accounting for symmetry
  if(!directed){
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind & upper.tri(s_ind,diag=TRUE)
    # linear mapping and projection basis
    Gd <- Sym_span(s_ind,s_ind_data)
    UV <- pracma::orth(crossprod(Gd,(Rproj %x% Lproj)[c(s_ind),]))
  }
  else{
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind
    # projection basis
    UV <- pracma::orth((Rproj %x% Lproj)[c(s_ind),])
  }
  # centering
  if(centered){
    UV <- center_orth(UV)
  }
  # final projection dimensions
  ss <- nrow(UV)
  dd <- ncol(UV)
  # projected data
  X <- crossprod(UV,sapply(A,function(x){x[s_ind_data]}))
  Y <- crossprod(UV,sapply(B,function(x){x[s_ind_data]}))
  # projected data means
  Xbar <- rowMeans(X)
  Ybar <- rowMeans(Y)
  # basic F test on augmented data assuming equal variance, independence
  # sse values for numerator
  sse_num <- sum((Xbar - Ybar)^2)
  df_num <- dd
  # adjustment for tilde, otherwise use projected data
  if(var_type=='orth'){
    Xa <- (diag(ss) - tcrossprod(UV)) %*% sapply(A,function(x){x[s_ind_data]})
    Ya <- (diag(ss) - tcrossprod(UV)) %*% sapply(B,function(x){x[s_ind_data]})
    sse_den <- sum((rowMeans(Xa) - rowMeans(Ya))^2)
    df_den <- ss - dd
  }
  else{
    sse_den <- (1/m1 + 1/m2)*(sum((X - Xbar)^2) + sum((Y - Ybar)^2))
    df_den <- dd*(m1+m2-2)
  }
  # statistic, dfs, pvalue
  out <- list()
  out$stat <- (sse_num/df_num)/(sse_den/(df_den))
  out$dfs <- df_num
  out$pval <- stats::pf(out$stat,df1=df_num,df2=df_den,lower.tail=FALSE)
  return(out)
}
