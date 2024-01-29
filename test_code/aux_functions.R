# auxilliary functions for data analysis

# expit function
expit <- function(x){
  1/(1+exp(-x))
}

# logistic regression variance function
vexpit <- function(x){
  expit(x)*(1-expit(x))
}

# logit function
logit <- function(x){
  log((x/(1-x)))
}

# replace diagonal with NAs
hollow <- function(M){
  temp <- M
  diag(temp) <- 0
  return(temp)
}

# replace diagonal with NAs
hollowna <- function(M){
  temp <- M
  diag(temp) <- NA
  return(temp)
}

# replace diagonal with 1s
hollowone <- function(M){
  temp <- M
  diag(temp) <- 1
  return(temp)
}

# basic exact F test for rectangles
# adjustment to ignore diagonal, account for symmetry, different number of
# networks in each group
basic_ftest <- function(A1,A2,hrow,hcol,alpha){
  # get small block
  # dimensions
  n <- dim(A1)[1]
  m1 <- dim(A1)[3]
  m2 <- dim(A2)[3]
  # find square and off-diagonal indices to keep
  # get the square part of the subset
  sqi <- intersect(hrow,hcol)
  sq <- matrix(FALSE,n,n)
  sq[sqi,sqi] <- TRUE
  sqvec <- c(as.logical(sq*upper.tri(sq,diag=FALSE)))
  # then the off-diagonal part
  od <- matrix(FALSE,n,n)
  od[hrow,hcol] <- TRUE
  od[sqi,sqi] <- FALSE
  odvec <- c(od)
  # number of observations for OLS
  N <- sum(c(sqvec,odvec))
  # terminate if either m is one
  if(m1==1 || m2==1){
    return(list(stat=NA,dfs=NA,pval=NA,rej=NA,conv=1))
  }
  # populate design matrix
  # populate full response matrices
  Y1 <- matrix(NA,m1,N)
  Y2 <- matrix(NA,m2,N)
  for(kk in 1:m1){
    Y1[kk,] <- c(A1[,,kk][sqvec],A1[,,kk][odvec])
  }
  for(kk in 1:m2){
    Y2[kk,] <- c(A2[,,kk][sqvec],A2[,,kk][odvec])
  }
  # sse values
  sse_null <- sum((scale(rbind(Y1,Y2),scale=F))^2)
  sse_alt <- sum((scale(Y1,scale=F))^2) + sum((scale(Y2,scale=F))^2)
  # statistic, dfs, pvalue
  stat <- ((sse_null - sse_alt)/N)/(sse_alt/(N*(m1+m2-2)))
  dfs <- N
  pval <- pf(stat,df1=dfs,df2=N*(m1+m2-2),lower.tail=FALSE)
  # determine rejection
  rej <- as.integer(pval < alpha)
  return(list(stat=stat,dfs=dfs,pval=pval,rej=rej))
}

# rank-truncated psuedoinverse
tpinv <- function(M,r){
  temp <- RSpectra::svds(M,r)
  dinv <- 1/temp$d
  out <- t(temp$u %*% (dinv * t(temp$v)))
  return(out)
}

# projection estimates from T
subspace_T <- function(A1bar,A2bar,d,hrow,hcol){
  # estimate blocks
  Chat <- A1bar[ hrow,-hcol] - A2bar[ hrow,-hcol]
  Rhat <- A1bar[-hrow, hcol] - A2bar[-hrow, hcol]
  Dhat <- A1bar[-hrow,-hcol] - A2bar[-hrow,-hcol]
  # estimate T
  That <- Chat %*% tpinv(Dhat,d) %*% Rhat
  # esimate subspaces
  temp <- svd(That)
  Lproj <- temp$u[,1:d]
  Rproj <- temp$v[,1:d]
  return(list(Lproj=Lproj,Rproj=Rproj))
}

# projection F test
# A1 and A2 are nxnxm arrays
# need to account for
# 1. diagonal
# 2. symmetry (follow discussion for non-rectangular hypothesis sets)
# 3. group-specific m's (reformulate statistic?)
proj_ftest_augment <- function(A1,A2,
                               hrow,hcol,alpha,
                               Lproj,Rproj,
                               on_diagonal=TRUE){
  # dimensions
  n <- dim(A1)[1]
  dL <- ncol(Lproj)
  dR <- ncol(Rproj)
  rr <- length(hrow)
  cc <- length(hcol)
  m1 <- dim(A1)[3]
  m2 <- dim(A2)[3]
  # test
  X <- matrix(NA,m1,dL*dR)
  Y <- matrix(NA,m2,dL*dR)
  for(kk in 1:m1){
    X[kk,] <- c(t(Lproj) %*% (A1[hrow,hcol,kk] %*% Rproj))
  }
  for(kk in 1:m2){
    Y[kk,] <- c(t(Lproj) %*% (A2[hrow,hcol,kk] %*% Rproj))
  }
  if(on_diagonal){
    X <- X[,c(upper.tri(diag(dL),diag=TRUE))]
    Y <- Y[,c(upper.tri(diag(dL),diag=TRUE))]
  }
  # remove columns which are repeated due to symmetry

  # basic F test on augmented data assuming equal variance, independence
  # sse values
  if((dL*dR)==1){
    dfs <- 1
    sse_null <-  sum((scale(c(X,Y),scale=F))^2)
  }
  else{
    dfs <- ncol(X)
    sse_null <- sum((scale(rbind(X,Y),scale=F))^2)
  }
  sse_alt <- sum((scale(X,scale=F))^2) + sum((scale(Y,scale=F))^2)
  # statistic, dfs, pvalue
  out <- list()
  out$stat <- (sse_null - sse_alt)/(sse_alt/(m1+m2-2))
  out$dfs <- dfs
  out$pval <- pf(out$stat,df1=dfs,df2=dfs*(m1+m2-2),lower.tail=FALSE)
  out$rej <- as.integer(out$pval < alpha)

  # instead try hotelling test with unequal variances
  # temp <- Hotelling::hotelling.test(X,Y,var.equal=FALSE)
  # out <- list()
  # out$stat <- temp$stats$statistic
  # out$dfs <- temp$stats$df[1]
  # out$pval <- temp$pval
  # out$rej <- as.integer(out$pval < alpha)
  return(out)
}

proj_ftest <- function(A1,A2,
                       hrow,hcol,alpha,
                       Lproj,Rproj,
                       tilde=FALSE){
  # dimensions
  n <- dim(A1)[1]
  dL <- ncol(Lproj)
  dR <- ncol(Rproj)
  rr <- length(hrow)
  cc <- length(hcol)
  m1 <- dim(A1)[3]
  m2 <- dim(A2)[3]
  # reset tilde if m is 1
  if(m1==1 || m2==1){
    tilde <- TRUE
  }
  # find square and off-diagonal indices to keep
  # get the square part of the subset
  sqi <- intersect(hrow,hcol)
  sq <- matrix(FALSE,n,n)
  sq[sqi,sqi] <- TRUE
  sqvec <- c(as.logical(sq*upper.tri(sq,diag=FALSE)))
  # then the off-diagonal part
  od <- matrix(FALSE,n,n)
  od[hrow,hcol] <- TRUE
  od[sqi,sqi] <- FALSE
  odvec <- c(od)
  # number of observations for OLS
  nn <- sum(c(sqvec,odvec))
  dd <- min(nn,dL*dR)
  # populate design matrix
  if(dd==nn){
    Z <- diag(nn)
  }
  else{
    Z <- matrix(NA,nn,dd)
    # dummy full projections
    Lproj_ext <- matrix(0,n,dL)
    Lproj_ext[hrow,] <- Lproj
    Rproj_ext <- matrix(0,n,dR)
    Rproj_ext[hcol,] <- Rproj
    cc <- 1
    for(rr in 1:dL){
      for(ss in 1:dR){
        temp <- Lproj_ext[,rr] %o% Rproj_ext[,ss]
        Z[,cc] <- c(temp[sqvec],temp[odvec])
        cc <- cc+1
      }
    }
  }
  # degrees of freedom
  nu1 <- dd
  if(tilde){
    nu2 <- (m1+m2)*(nn - dd)
  }
  else{
    nu2 <- (m1+m2-2)*(nn - dd)
  }
  # data matrices
  # sample 1
  Y1 <- matrix(NA,nn,m1)
  for(kk in 1:m1){
    Y1[,kk] <- c(A1[,,kk][sqvec],A1[,,kk][odvec])
  }
  # sample2
  Y2 <- matrix(NA,nn,m2)
  for(kk in 1:m2){
    Y2[,kk] <- c(A2[,,kk][sqvec],A2[,,kk][odvec])
  }
  # kronecker o/n basis (nn x dd) # needs to be orthonormal
  UV <- svd(Z)$u
  # difference subspace, weighted by relative sample sizes
  W <- rbind(UV*sqrt(m2/(m1+m2)),-UV*sqrt(m1/(m1+m2))) # 2nn x dd, still o/n
  # complement projection operator (nn x dd), idempotent
  N <- diag(nn) - tcrossprod(UV)
  # centering matrix
  one_m1 <- rep(1,m1)
  one_m2 <- rep(1,m2)
  if(!tilde){
    center_m1 <- diag(m1) - tcrossprod(one_m1)/m1
    center_m2 <- diag(m2) - tcrossprod(one_m2)/m2
    #center_M <- diag(m1+m2) - tcrossprod(rep(1,m1+m2))/(m1+m2)
  }
  # calculate numerator
  numvec <- crossprod(W,c((Y1 %*% one_m1)/sqrt(m1), (Y2 %*% one_m2)/sqrt(m2)))
  num <- nu2*sum(numvec^2)
  # calculate denominator
  if(tilde){
    denmat <- crossprod(N,cbind(Y1,Y2))
  }
  else{
    denmat <- crossprod(N,cbind(Y1 %*% center_m1,Y2 %*% center_m2))
    #denmat <- crossprod(N,cbind(Y1,Y2) %*% center_M)
  }
  den <- nu1*sum(denmat^2)
  # compile results and return
  out <- list()
  out$stat <- num/den
  out$dfs <- nu1
  out$pval <- pf(out$stat,nu1,nu2,lower.tail=FALSE)
  out$rej <- as.integer(out$pval < alpha)
  return(out)
}

# basic exact with binomial proportion testing
basic_binprop <- function(A1,A2,hrow,hcol,
                                  alpha){
  # dimensions
  n <- dim(A1)[1]
  m1 <- dim(A1)[3]
  m2 <- dim(A2)[3]
  # terminate if m=1
  if(m1==1 || m2==1){
    return(list(stat=0,dfs=N,pval=1,rej=0))
  }
  # find square and off-diagonal indices to keep
  # get the square part of the subset
  sqi <- intersect(hrow,hcol)
  sq <- matrix(FALSE,n,n)
  sq[sqi,sqi] <- TRUE
  sqvec <- c(as.logical(sq*upper.tri(sq,diag=TRUE))) # keep diagonal for cell data
  # then the off-diagonal part
  od <- matrix(FALSE,n,n)
  od[hrow,hcol] <- TRUE
  od[sqi,sqi] <- FALSE
  odvec <- c(od)
  # number of observations for OLS
  N <- sum(c(sqvec,odvec))
  # populate design matrix
  # populate response vectors
  y1 <- y2 <- rep(0,N)
  for(kk in 1:m1){
    y1 <- y1 + c(A1[,,kk][sqvec],A1[,,kk][odvec])
  }
  for(kk in 1:m2){
    y2 <- y2 + c(A2[,,kk][sqvec],A2[,,kk][odvec])
  }
  # combined chi square stat
  stat <- 0
  dfs <- 0
  for(ii in 1:N){
    suppressWarnings(test <- prop.test(x=c(y1[ii],y2[ii]),
                                       n=c(m1,m1),
                                       correct=FALSE))
    chi <- as.vector(test$stat)
    # update statistic
    if(!is.nan(chi)){
      stat <- stat + chi
    }
    # update dfs
    dfs <- dfs + as.vector(test$parameter)
  }
  # dfs, pvalue
  pval <- pchisq(stat,df=dfs,lower.tail=FALSE)
  # determine rejection
  rej <- as.integer(pval < alpha)
  return(list(stat=stat,dfs=dfs,pval=pval,rej=rej))
}

# projection W test
# A1 and A2 are nxnxm arrays
proj_wtest <- function(A1,A2,
                       hrow,hcol,
                       Lproj,Rproj,
                       tilde=FALSE,
                       on_diagonal=FALSE,
                       add_intercept=FALSE,
                       misspec=TRUE){
  # dimensions
  rr <- length(hrow)
  cc <- length(hcol)
  dL <- ncol(Lproj)
  dR <- ncol(Rproj)
  m1 <- dim(A1)[3]
  m2 <- dim(A2)[3]
  # reset tilde if m is 1
  if(m1==1 || m2==1){
    tilde <- FALSE
  }
  # manually construct diagonal basis or use current code
  if(on_diagonal){
    sqvec <- c(upper.tri(diag(rr),diag=TRUE))
    UV <- NULL
    for(s1 in 1:dL){
      for(s2 in 1:s1){
        if(s1==s2){
          UV <- cbind(UV,c(Lproj[,s1] %o% Rproj[,s2])[sqvec])
        }
        else{
          UV <- cbind(UV,c(c(Lproj[,s1] %o% Rproj[,s2]) + c(Lproj[,s1] %o% Rproj[,s2]))[sqvec])
        }
      }
    }
    if(add_intercept){
      UV <- cbind(1/sqrt(sum(sqvec)),scale(UV,scale=F))
    }
    dd <- ncol(UV)
  }
  else{
    # kronecker o/n basis (rc x (2d)^2)
    UV <- Rproj %x% Lproj
    if(add_intercept){
      UV <- cbind(1/sqrt(sum(rr*cc)),scale(UV,scale=F))
    }
    dd <- ncol(UV)
  }
  # common subspace
  X <- rbind(UV,UV)
  # difference subspace
  W <- rbind(UV,-UV)
  # final design
  Xall <- cbind(X,W)
  # populate full successes
  if(on_diagonal){
    s_vec <- rep(0,2*sum(sqvec))
    for(kk in 1:m1){
      s_vec <- s_vec + c(A1[hrow,hcol,kk][sqvec],rep(0,sum(sqvec)))
    }
    for(kk in 1:m2){
      s_vec <- s_vec + c(rep(0,sum(sqvec)),A2[hrow,hcol,kk][sqvec])
    }
  }
  else{
    s_vec <- rep(0,2*rr*cc)
    for(kk in 1:m1){
      s_vec <- s_vec + c(A1[hrow,hcol,kk],rep(0,rr*cc))
    }
    for(kk in 1:m2){
      s_vec <- s_vec + c(rep(0,rr*cc),A2[hrow,hcol,kk])
    }
  }
  # calculate failures
  f_vec <- c(rep(m1,length(s_vec)/2),rep(m2,length(s_vec)/2)) - s_vec
  # fit logistic regression
  if(tilde){
    moda <- glm(cbind(s_vec,f_vec)~Xall-1,family=quasibinomial())
  }
  else{
    moda <- glm(cbind(s_vec,f_vec)~Xall-1,family=binomial())
  }
  # get logistic regression coefficients and linear predictor
  if(add_intercept){
    gamma1 <- moda$coefficients[1:dd]
    gamma2 <- moda$coefficients[-(1:(dd+1))]
    #print(gamma2)
    lmhat <- c(UV %*% gamma1)
  }
  else{
    gamma1 <- moda$coefficients[1:dd]
    gamma2 <- moda$coefficients[-(1:dd)]
    #print(gamma2)
    lmhat <- c(UV %*% gamma1)
  }
  # estimate covariance
  # saturated model predictor
  conserve <- 0.5
  prop_vec <- (s_vec[1:(length(s_vec)/2)] + s_vec[-(1:(length(s_vec)/2))])/(m1+m2)
  prop_vec_reg <- sign(prop_vec - 0.5)*pmax(abs(prop_vec - 0.5) - (conserve/(m1+m2)),0) + 0.5
  smhat <- logit(prop_vec_reg)
  # covariances
  if(add_intercept){
    Ghat <- 2*(t(matrix(UV[,-1],ncol=dd-1)) %*% (diag(vexpit(lmhat)) %*% matrix(UV[,-1],ncol=dd-1)))
    Fhat <- 2*(t(matrix(UV[,-1],ncol=dd-1)) %*% (diag(vexpit(smhat)) %*% matrix(UV[,-1],ncol=dd-1)))
  }
  else{
    Ghat <- 2*(t(UV) %*% (diag(vexpit(lmhat)) %*% UV))
    Fhat <- 2*(t(UV) %*% (diag(vexpit(smhat)) %*% UV))
  }
  # full inverse covariance
  if(misspec){
    iVhat <- Ghat %*% (solve(Fhat) %*% Ghat)
  }
  else{
    iVhat <- Ghat
  }
  # stat numerator
  num <- ((m1+m2)/2)*sum(gamma2 * (iVhat %*% gamma2))
  # overdispersion estimate as denominator
  den <- summary(moda)$dispersion
  # compile results and return
  out <- list()
  out$stat <- num/den
  out$dfs <- dd
  out$pval <- pchisq(out$stat,dd,lower.tail=FALSE)
  return(out)
}

# scratch code: constructing symmetry constraints

# # projections
# U <- svd(matrix(rnorm(20),ncol=2))$u
# V <- svd(matrix(rnorm(20),ncol=2))$u
#
# # symmetry constraints on first block
# G <- matrix(0,200,45)
# kk <- 1
# for(ii in 1:10){
#   for(jj in 1:ii){
#     if(jj < ii){
#       gtemp <- matrix(0,10,10)
#       gtemp[ii,jj] <- 1
#       gtemp[jj,ii] <- -1
#       G[,kk] <- c(gtemp,rep(0,100))
#       kk <- kk+1
#     }
#   }
# }
#
# # convert to a constraint on the coefficients
# Gtil <- crossprod(G,rbind(U,V) %x% U)
# temp <- svd(Gtil)
# Gnull <- temp$v[,temp$d < 1e-6]




