library(Matrix)
library(irlba)
library(RMTstat)

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
#' @param A cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#' @param B cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#' @param sig significance level for acceptance of null hypothesis
#' @param r r communities in G
#'
#' @return acceptance/rejection decision & p-value for a normal distribution
#' @export
#'
#' @examples
#' simu = twosamp_twoblock(n = 4, m = 4, p = 0.3, q = 0.3, epsilon = 0.7)
#' test_result = Asymp_TW(simu$A_G, simu$A_H, 0.05, 100)
#'
Asymp_TW <- function(A, B, sig, r) {
  n = dim(A[[1]])[1]

  if (length(r)==1){
    idx = spectral_clus(A, B,r);
  }
  else{
    idx = r
    r = max(idx) #k
  }

  #compute C
  C = A[[1]] - B[[1]]

  for (i in 1:r){
    for (j in 1:r) {
      if (i != j){
        curr = A[[1]][idx == i, idx == j]
        Pij = mean(curr)

        curr = B[[1]][idx == i, idx == j]
        Qij = mean(curr)
      }
      else{
        curr = A[[1]][idx == i, idx == j]
        Pij = sum(curr) / (dim(curr)[1] * (dim(curr)[1] - 1)) # Pij might = x / 0

        curr = B[[1]][idx == i, idx == j]
        Qij = sum(curr) / (dim(curr)[1] * (dim(curr)[1] - 1))
      }

      numer = C[idx == i, idx == j]
      denom = sqrt((n - 1) * (Pij) * (1 - Pij) + Qij * (1 - Qij))

      if (identical(denom, numeric(0))){
        denom = 1e-5
      }
      if (denom == 0){
        denom = 1e-5
      }
      C[idx == i, idx == j] = numer / denom
    }
  }
  # NA check
  C[is.na(C)] = 0
  #test statistics
  spectral.norm = irlba(C, 1)$d # the largest singular value
  test.stat = n ^ (2 / 3) * (spectral.norm - 2)
  p.val = ptw(test.stat, beta=1, lower.tail = FALSE, log.p = FALSE)
  test <- ifelse(p.val <= sig, 1, 0)

  return(c(test, p.val))
}


# Auxiliary functions

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
#' @param A a sparse EPAacency matrix
#' @param B a sparse EPAacency matrix
#' @param r r communities in G
#'
#' @return returns a vector of length n containing the cluster assignment of each node.
#' @export
#'
#' @examples
#' spectral_clus(simu$A_G, simu$A_H, r = 3)
spectral_clus = function(A, B, r){
  C = (A[[1]] + B[[1]])/2
  d = rowSums(C)
  d[d==0] = 1
  d = 1/sqrt(d)
  C = C * (d %*% t(d))
  vec = svd(C)$v[, 1:r]
  normv = sqrt(rowSums(vec^2))
  normv[normv==0] = 1

  a = matrix(normv, nrow=1, ncol=r, byrow=TRUE)
  output = vec

  for(i in 1:nrow(vec)){
    output[i, ] = vec[i, ] / a
  }
  idx = kmeans(output, centers=r)$cluster
  return(idx)
}

