#### helper functions ####

# laplacian spectrum computation for one adjacency matrix
# Compute the eigenvalue spectrum of either the unnormalized Laplacian
# (\eqn{L = D - A}) or the normalized Laplacian
# (\eqn{L_{\text{sym}} = I - D^{-1/2} A D^{-1/2}}).
# Supports truncation to the first/last \code{k} eigenvalues and optional L2 normalization.
# keep: Number of eigenvalues to keep (default = NULL = all)

# which argument is deprecated, always takes the smallest as in SMou_2025/Simulation_code
# normalization of spectrum is deprecated, always normalizes as in Huang et al.
laplacian_spectrum <- function(A,laplacian = 'normalized',keep = NULL) {
  # initial cleaning/option matching
  laplacian <- match.arg(laplacian,c("normalized", "unnormalized"))
  # compute normalized or unnormalized laplacian
  if (laplacian == "normalized") {
    deg <- rowSums(A)
    inv_sqrt_deg <- 1 / sqrt(pmax(deg, 1e-12))
    Dm12 <- diag(inv_sqrt_deg)
    L <- diag(nrow(A)) - Dm12 %*% A %*% Dm12
  } else {
    D <- diag(rowSums(A))
    L <- D - A
  }
  # compute spectrum, sort decreasing
  ev <- eigen(L, only.values = TRUE, symmetric = TRUE)$values
  ev <- sort(ev, decreasing = FALSE)
  # truncate spectrum if keep < n
  if (!is.null(keep)) {
    k <- min(keep, length(ev))
    ev <- ev[seq_len(k)]
  }
  return(ev)
}

# computation of context vectors from a list of spectra, either svd or mean aggregation
context_vector <- function(spectra_list, method = 'svd') {
  method <- match.arg(method,c('svd','mean'))
  mat <- do.call(cbind, spectra_list)
  if (method == 'svd') {
    return(c(irlba::irlba(mat,1)$u))
  } else {
    return(rowMeans(mat))
  }
}

#### main function ####
# NOTE: distance argument is deprecated, always uses default cos_dis as in Huang et al.

#' Graph Laplacian Spectral Anomaly Scores for Changepoint Detection
#'
#' Computes 2 normalized anomaly scores based on Laplacian spectra,
#' described in \href{https://doi.org/10.1145/3394486.3403077}{Huang et al., (2020)}.
#' for an ordered sequence of undirected networks with no self-loops.
#' 1. \code{Z}: raw anomaly score (cosine distance from short-term or long-term context vectors).
#' 2. \code{Zstar}: change in raw anomaly score.
#'
#' @param A A list of matrices, ordered sequence of adjacency matrices.
#' @param context_method Either \code{'svd'} or \code{'mean'}, method for constructing spectral context vectors. Defaults to \code{'svd'}.
#' @param context_short Integer, short-term context window size, defaults to \code{5}.
#' @param context_long Integer, long-term context window size, defaults to \code{25}.
#' @param laplacian Either \code{'normalized'} or \code{'unnormalized'}, which graph Laplacian to use. Defaults to \code{'normalized'}.
#' @param keep Integer, how many eigenvalues to keep from Laplacian spectrum.
#'
#' @return A \code{length(A)}\eqn{\times}\code{2} matrix with named columns, containing the 2 anomaly scores.
#'
#' @export
#'
#' @examples
#' CP <- list(type='kidneyegg',t_start=30,n_egg=20,p=0.5,q=0.8)
#' A <- Simulate_CP(40,40,CP=CP)$A
#' test <- CP_laplacian(A)
CP_laplacian <- function(A,context_method = 'svd',
                         context_short=5,context_long=25,
                         laplacian = 'normalized',
                         keep = NULL) {
  # initial argument checking
  context_method <- match.arg(context_method,c("svd", "mean"))
  laplacian <- match.arg(laplacian,c("normalized", "unnormalized"))
  # dimensions and data cleaning
  A <- checklist(A)
  m <- length(A)
  Z <- Zstar <- rep(NA,m)
  # compute spectra
  spectra <- lapply(A,laplacian_spectrum,laplacian=laplacian,keep=keep)
  # long window context vector
  long_vec <- context_vector(spectra[seq_len(context_long)], method = context_method)
  # starting time for computation
  start_t <- context_long + 1
  # populate Z vector
  if (start_t <= m) {
    for (tt in start_t:m) {
      short_vec <- context_vector(spectra[(tt - context_short):(tt - 1)], method = context_method)
      # compute Z-scores
      Zs <- cos_dis(spectra[[tt]], short_vec)
      Zl <- cos_dis(spectra[[tt]], long_vec)
      Z[tt] <- max(Zs, Zl)
    }
  }
  else{
    stop('Context window too long to compute scores')
  }
  # populate Zstar vector
  idx <- which(!is.na(Z))
  if (length(idx) >= 2) {
    Zstar[idx[-1]] <- pmax(Z[idx[-1]] - Z[idx[-length(idx)]], 0)
  }
  # compile matrix and return
  out <- cbind(Z,Zstar)
  colnames(out) = c('Z','Zstar'); rownames(out) <- 1:m
  # trim NA entries
  out_trim <- trimNA(out)
  return(out_trim)
}


