# Utility functions for nettest

#### R utilities ####

# helper: hollowize a square matrix (set diagonal entries to zero)
hollowize <- function(A){
  A - diag(diag(A))
}

# helper: convert array to list of aligned matrices
array_to_list <- function(A,self_loops=TRUE){
  m <- dim(A)[3]
  if(self_loops){
    lapply(1:m,function(kk){A[,,kk]})
  }
  else{
    lapply(1:m,function(kk){hollowize(A[,,kk])})
  }
}

# helper: convert list of aligned matrices to 3d array
list_to_array <- function(B){
  n <- nrow(B[[1]])
  m <- length(B)
  array(unlist(B),c(n,n,m))
}

# probability clipping function
pclip <- function(p,eps){
  pmax(pmin(p,1-eps),eps)
}

# convert a categorial vector C to a binary membership matrix Z
C_to_Z <- function(C,K){
  n <- length(C)
  Z <- matrix(0,n,K)
  Z[cbind(seq_len(n),C)] <- 1
  return(Z)
}

# convert a binary membership matrix Z to a categorical vector C
Z_to_C <- function(Z){
  apply(Z,1,function(v){which.max(v)})
}

#### spectral clustering ####

# Spectral Clustering used in TW test -------------------------------------

# Example: # clusters <- spectral_clus(simu$A_G, simu$A_H, r = 3)

#' spectral_clus
#'
#'
#' This function computes the cluster assignment of each node, serving as an auxiliary functions for the function: Asymp_TW
#'
#'
#'
#'
#' @param C a matrix
#' @param r r communities in G
#'
#' @return returns a vector of length n containing the cluster assignment of each node.
#' @export
#'
#' @examples
#' C = tcrossprod(1:10)
#' spectral_clus(C, r = 3)
spectral_clus = function(C, r){
  d = rowSums(C)
  d[d==0] = 1
  d = 1/sqrt(d)
  C = C * (d %*% t(d))
  vec = svd(C)$v[, 1:r]
  normv = sqrt(rowSums(vec^2))
  normv[normv==0] = 1

  a = matrix(rep(normv,r), ncol=r)
  output = vec/a

  idx = stats::kmeans(output, centers=r)$cluster
  return(idx)
}

# populating a matrix with its block averages
block_avg <- function(M,c){
  Mbar <- matrix(0,nrow(M),ncol(M))
  r <- max(c)
  apply(expand.grid(1:r,1:r),1,function(s){
    ix <- which(c==s[1])
    iy <- which(c==s[2])
    Mbar[ix,iy] <<- mean(M[ix,iy])
  })
  return(Mbar)
}

#### mesoscale testing ####

# logit function
logit <- function(x){
  log((x/(1-x)))
}

# expit function
expit <- function(x){
  1/(1+exp(-x))
}

# logistic regression variance function
vexpit <- function(x){
  expit(x)*(1-expit(x))
}

# rank-truncated psuedoinverse
Tpinv <- function(M,r){
  temp <- irlba::irlba(M,r)
  dinv <- 1/temp$d
  out <- t(temp$u %*% (dinv * t(temp$v)))
  return(out)
}

# calculate the pseudoinverse matrix for enforcing symmetry constraints
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
  if(d==1){
    return(list(Lproj=matrix(Lproj,ncol=1),Rproj=matrix(Rproj,ncol=1)))
  }
  else{
    return(list(Lproj=Lproj,Rproj=Rproj))
  }
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

#### OMNI test ####

# genIERGraph <- function(m, model) {
#   A <- list()
#   for (i in 1:m) {
#     EA <- upper.tri(model$P)
#     A1 <- matrix(runif(length(EA)) < EA, nrow = ncol(model$P))
#     A1 <- A1 + t(A1)  # Ensure symmetry by adding the transpose
#     A[[i]] <- A1
#   }
#   return(A)
# }

get_omnibus_matrix_sparse <- function(matrices) {
  rows <- list()
  # Iterate over each column
  for (column_index in seq_along(matrices)) {
    current_row <- list()

    for (row_index in seq_along(matrices)) {
      if (row_index == column_index) {
        # we are on the diagonal, we do not need to perform any calculation and instead add the current matrix
        # to the current_row
        current_row[[row_index]] <- matrices[[column_index]]
      } else {
        # otherwise we are not on the diagonal and we average the current_matrix with the matrix at row_index
        # and add that to our current_row
        matrices_averaged <- (matrices[[column_index]] + matrices[[row_index]]) * 0.5
        current_row[[row_index]] <- matrices_averaged
      }
    }

    result <- matrix(0, nrow=0, ncol=ncol(current_row[[1]]))

    # Iterate through the current_row list and vertically concatenate matrices into the result matrix
    for (i in 1:length(current_row)) {
      result <- rbind(result, current_row[[i]])
    }

    # row
    rows[[column_index]] <- result
  }

  # Combine rows to create the omnibus matrix
  omnibus_matrix <- matrix(0, nrow = nrow(rows[[1]]), ncol=0)

  # Iterate through the rows list and horizontally concatenate matrices into the result matrix
  for (i in 1:length(rows)) {
    omnibus_matrix <- cbind(omnibus_matrix, rows[[i]])
  }

  return(omnibus_matrix)
}

# NOTE: may use this later
#' Extract Community Labels from Block Membership Matrix
#'
#' @description
#' Converts a binary block membership matrix \eqn{Z} into a vector of community
#' labels. Each row of \eqn{Z} is assumed to be a one-hot encoding of a node’s
#' community membership (i.e., exactly one entry per row is equal to 1).
#'
#' @param Z A binary membership matrix of size \eqn{n \times K}, where
#'   \eqn{n} is the number of nodes and \eqn{K} is the number of communities.
#'
#' @return An integer vector of length \eqn{n} giving the community label
#'   (from 1 to \eqn{K}) for each node.
#'
#' @examples
#' # Example: 4 nodes, 2 communities
#' Z <- matrix(c(1,0,
#'               1,0,
#'               0,1,
#'               0,1), nrow = 4, byrow = TRUE)
#' get_labels_from_Z(Z)
#' # Returns: c(1, 1, 2, 2)
#'
# get_labels_from_Z <- function(Z) {
#   apply(Z, 1, function(row) which(row == 1))
# }

