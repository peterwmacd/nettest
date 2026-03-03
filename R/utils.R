# Utility functions for nettest

#### exported utilities ####

#' (Normalized) Spectral Clustering
#'
#' This function computes the cluster assignment based on the normalized singular vectors of a matrix \code{M}.
#'
#' @param M a square \code{n}\eqn{\times}\code{n} matrix
#' @param d number of clusters
#'
#' @return returns a vector of length \code{n} containing the integer-valued cluster assignments.
#'
#' @export
#'
#' @examples
#' M = tcrossprod(1:10)
#' spectral_clust(M, d = 3)
spectral_clust = function(M,d){
  # normalize input matrix
  theta = rowSums(M)
  theta[theta==0] = 1
  theta = 1/sqrt(theta)
  Mnorm = M * (theta %*% t(theta))
  # get normalized top d singular vectors
  svec = irlba::irlba(Mnorm,d)$v
  normv = sqrt(rowSums(svec^2))
  normv[normv==0] = 1
  svec_norm = svec/matrix(rep(normv,d), ncol=d)
  # apply k-means clustering
  C = stats::kmeans(svec_norm, centers=d)$cluster
  return(C)
}

#### other utilities ####

# convert array to list of aligned matrices
array_to_list <- function(A,self_loops=TRUE){
  m <- dim(A)[3]
  if(self_loops){
    lapply(1:m,function(kk){A[,,kk]})
  }
  else{
    lapply(1:m,function(kk){hollowize(A[,,kk])})
  }
}

# average over blocks of a matrix (row-wise) w/kronecker
# for combined omni embedding
avg_block <- function(M,nblock){
  # dimension
  n <- nrow(M)/nblock
  # initialize
  crossprod((rep(1,nblock) %x% diag(n)),M) / nblock
}

# merging/averaging the first two communities from a B-matrix
B_merge <- function(B){
  # dimension
  K <- ncol(B)
  # linear operator to merge
  M <- rbind(c(0.5,rep(0,K-2)),diag(K-1)); M[2,1] <- 0.5
  # apply transformation
  crossprod(M , (B %*% M))
}

# populating a matrix with its block averages
# full=TRUE produces an nxn matrix, full=FALSE produces a KxK matrix
block_avg <- function(M,Z,full=TRUE){
  if(full){
    D <- Z %*% solve(crossprod(Z),t(Z))
  }
  else{
    D <- solve(crossprod(Z),t(Z))
  }
  D %*% tcrossprod(M,D)
}

# Bayesian estimator of mean for proportion data
# equivalent to posterior mean under a Unif(0,1) prior
bmean <- function(v){
  (sum(v) + 1)/(length(v) + 2)
}

# column-wise comparison to a bootstrap null distribution, with a continuity correction
# for one p-value
bootp <- function(t,t0){
  # dimension
  n0 <- length(t0)
  # bootstrap p
  (sum(t0 >= t) + 0.5) / n0
}

# for multiple p-values
bootp_multi <- function(tvec,t0mat){
  # dimension
  q <- ncol(t0mat)
  # loop over bootstrap p
  pvec <- numeric(q)
  for(jj in 1:q){
    pvec[jj] <- bootp(tvec[jj],t0mat[,jj])
  }
  return(pvec)
}

# convert a categorial vector C to a binary membership matrix Z
C_to_Z <- function(C,K){
  n <- length(C)
  Z <- matrix(0,n,K)
  Z[cbind(seq_len(n),C)] <- 1
  return(Z)
}

# 'centering' function for an n times d orthonormal matrix, ie make it
# strictly orthogonal to the 1's vector
center_orth <- function(W){
  n <- nrow(W)
  cmat <- diag(n) - matrix(1,n,n)/n
  Wc <- cmat %*% W
  Wc_orth <- pracma::orth(Wc)
  return(Wc_orth)
}

# check that an object is a list, if not make it a length 1 list
checklist <- function(L){
  if(is.list(L)){
    return(L)
  }
  else{
    return(list(L))
  }
}

# cosine angular distance between two vectors (performs normalization internally)
cos_dis <- function(a, b) {
  denom <- sqrt(sum(a*a)) * sqrt(sum(b*b))
  if(denom == 0){
    return(0)
  }
  else{
    return(1 - abs(sum(a * b) / denom))
  }
}

# compute pairwise euclidean distances between the rows of two matrices
eucdist <- function(X,Y){
  # dimensions
  n <- nrow(X)
  # distance matrix
  as.matrix(stats::dist(rbind(X,Y)))[1:n,-(1:n)]
}

# expit function
expit <- function(x){
  1/(1+exp(-x))
}

# hollowize a square matrix (set diagonal entries to zero)
hollowize <- function(A){
  A - diag(diag(A))
}

# convert list of aligned matrices to 3d array
list_to_array <- function(B){
  n <- nrow(B[[1]])
  m <- length(B)
  array(unlist(B),c(n,n,m))
}

# convert list of aligned matrices to vectorized data
list_to_vecdata <- function(B){
  n <- nrow(B[[1]])
  m <- length(B)
  t(matrix(unlist(lapply(B,pracma::triu)),n*n,m))
}

# logit function
logit <- function(x){
  log((x/(1-x)))
}

# probability clipping function
pclip <- function(p,eps){
  pmax(pmin(p,1-eps),eps)
}

# elementwise perturbation of a matrix by Gaussian noise
perturb_mat <- function(M,sd){
  # perturb
  M + matrix(stats::rnorm(prod(dim(M)),sd=sd),nrow(M))
}

# extending Pi (SBM probability) vector by splitting the probabilities for the first community
Pi_split <- function(Pi){
  c(Pi[1]/2,Pi[1]/2,Pi[-1])
}

# procrustes alignment applied to Y (first argument) to align with X (second argument)
proc_align <- function(Y,X){
  # factorize crossproduct
  temp <- svd(crossprod(Y,X))
  # get transformation
  orth <- tcrossprod(temp$u,temp$v)
  Y_orth <- Y %*% orth
  return(Y_orth)
}

# symmetric Gaussian noise with dispersion (edge variance)
sym_noise <- function(n,dispersion=1){
  nhalf <- choose(n+1,2)
  E <- matrix(0,n,n)
  E[upper.tri(E,diag=TRUE)] <- stats::rnorm(nhalf,sd=sqrt(dispersion))
  E[lower.tri(E,diag=FALSE)] <- t(E)[lower.tri(E,diag=FALSE)]
  return(E)
}

# calculate the pseudoinverse matrix for enforcing symmetry constraints on a set
# of matrix indices. s_ind, s_ind_tri are matrices of booleans, describing the active
# index sets (overall and upper triangular resp.)
Sym_span <- function(s_ind,s_ind_tri){
  # reordered hypothesis set and upper triangle
  ind <- which(s_ind,arr.ind=TRUE)
  ind_tri <- which(s_ind_tri,arr.ind=TRUE)
  # G dimensions
  n_ind <- nrow(ind)
  n_ind_tri <- nrow(ind_tri)
  # construct G matrix
  G <- matrix(0,n_ind,n_ind_tri)
  for(kk in 1:n_ind){
    entry <- ind[kk,]
    tri_index <- which(((entry[1] == ind_tri[,1]) & (entry[2] == ind_tri[,2])) | ((entry[1] == ind_tri[,2]) & (entry[2] == ind_tri[,1])))
    G[kk,tri_index] <- 1
  }
  # (left) pseudoinverse
  Gd <- t(t(G)/colSums(G))
  return(Gd)
}

# rank-truncated psuedoinverse
Tpinv <- function(M,r){
  temp <- irlba::irlba(M,r)
  dinv <- 1/temp$d
  out <- t(temp$u %*% (dinv * t(temp$v)))
  return(out)
}

# remove rows of a matrix that are all NA
trimNA <- function(M){
  M[apply(M,1,function(x){!all(is.na(x))}),]
}

# populating a vector with its block averages
vblock_avg <- function(v,Z){
  Z %*% solve(crossprod(Z),crossprod(Z,v))
}

# populating a vector with its block standard deviations
vblock_sd <- function(v,C){
  # allocate space
  K <- max(C)
  w <- numeric(K)
  for(kk in seq_len(K)){
    ind <- (C==kk)
    w[kk] <- stats::sd(v[ind])
  }
  return(w)
}

# logistic regression variance function
vexpit <- function(x){
  expit(x)*(1-expit(x))
}

# merging the first two communities from a Z-matrix
Z_merge <- function(Z){
  # dimension
  K <- ncol(Z)
  # linear operator to merge
  M <- rbind(c(1,rep(0,K-2)),diag(K-1))
  # apply transformation
  Z %*% M
}

# convert a binary membership matrix Z to a categorical vector C
Z_to_C <- function(Z){
  apply(Z,1,function(v){which.max(v)})
}
