#### helper functions ####

# names for vector of dcsbm statistics
dcsbm_names <- function(K){
  namevec <- numeric(1 + choose(K+1,2) + K)
  namevec[1] <- 'pimax'
  kk <- 2
  for(ii in 1:K){
    for(jj in 1:ii){
      namevec[kk] <- paste0('B',jj,ii)
      kk <- kk+1
    }
  }
  for(ii in 1:K){
    namevec[kk] <- paste0('theta',ii)
    kk <- kk+1
  }
  return(namevec)
}

# compute DCSBM statisitcs based a adjacency matrix A, and community membership
# vector C
dcsbm_stats <- function(A,C){
  # dimensions
  n <- nrow(A)
  K <- max(C)
  # allocate space
  nstat <- 1 + choose(K+1,2) + K
  statvec <- numeric(nstat)
  # pimax from normalized spectral clustering
  Chat <- spectral_clust(A,K)
  statvec[1] <- max(table(Chat)/n)

  # local edge densities
  Z <- C_to_Z(C,K)
  B <- block_avg(A,Z,full=FALSE)
  statvec[2:(1 + choose(K+1,2))] <- B[upper.tri(B,diag=TRUE)]

  # local standard deviations of theta
  dvec <- rowSums(A)
  dmeans <- vblock_avg(dvec,Z)
  statvec[(2 + choose(K+1,2)):nstat] <- vblock_sd(dvec/dmeans,C)

  # return result
  return(statvec)
}

#### main function ####

#' DCSBM Statistics for Changepoint Detection
#'
#' Computes 9 graph summary statistics (or graph invariants),
#' described in \href{https://doi.org/10.1002/qre.2520}{Wilson et al., (2019)}.
#' for an ordered sequence of undirected networks with no self-loops.
#' 1. \code{pimax}: proportion of nodes in largest community, according to spectral clustering.
#' 2. \code{Bkl}: local densities for \eqn{1 \leq l \leq k \leq K}.
#' 3. \code{thetak}: sample standard deviations of block \eqn{k} degree-correction parameters.
#'
#' @param A A list of matrices, ordered sequence of adjacency matrices.
#' @param C A vector, categorical/integer community memberships for computing local statistics. Defaults to \code{NULL}.
#' @param K An integer, number of communities for normalized spectral clustering (based on the first adjacency matrix in the sequence). Defaults to \code{2}. Only used if \code{Z} is unspecified, otherwise \code{K = max(C)}.
#'
#' @return A \code{length(A)}\eqn{\times}\code{(1 + choose(K+1,2) + K)} matrix with named columns, with the summary statistics computed for each adjacency matrix.
#'
#' @export
#'
#' @examples
#' CP <- list(type='kidneyegg',t_start=5,n_egg=20,p=0.5,q=0.8)
#' data <- Simulate_CP(100,10,CP=CP)
#' A <- data$A; C <- Z_to_C(data$Z1)
#' test <- CP_dcsbm(A,C)
CP_dcsbm <- function(A,C=NULL,K=2){
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  m <- length(A)
  # populate C if unspecified (spectral clustering on initial snapshot), otherwise
  # override K
  if(is.null(C)){
    C <- spectral_clust(A[[1]],K)
  }
  else{
    K <- max(C)
  }
  # compute summaries
  summaries <- t(apply(list_to_array(A),3,dcsbm_stats,C=C))
  colnames(summaries) <- dcsbm_names(K)
  rownames(summaries) <- 1:m
  return(summaries)
}
