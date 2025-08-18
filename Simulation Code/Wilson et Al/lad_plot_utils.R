#' Compute Normalized Laplacian Singular Spectrum
#'
#' @param A Adjacency matrix (numeric, square)
#' @param k Optional: number of top singular values to retain
#'
#' @return A unit-length numeric vector of singular values
#' @export
get_laplacian_spectrum <- function(A, k = NULL) {
  L <- diag(rowSums(A)) - A
  sv <- svd(L, nu = 0, nv = 0)$d
  if (!is.null(k)) sv <- sv[1:k]
  return(sv / sqrt(sum(sv^2)))  # L2 normalize
}

#' Compute Context Vector from Spectra
#'
#' @param spectra_list List of normalized singular spectra (each a numeric vector)
#'
#' @return The first left singular vector of the stacked spectra matrix
#' @export
get_context_vector <- function(spectra_list) {
  mat <- do.call(cbind, spectra_list)
  svd(mat)$u[, 1]
}

#' Compute LAD Z and Z* Scores
#'
#' @param spectra List of normalized Laplacian spectra (from each time step)
#' @param s Size of the short-term window (default: 5)
#' @param l Size of the long-term window (default: 10)
#'
#' @return A list with components:
#'   \item{Z}{Raw anomaly scores (cosine dissimilarity)}
#'   \item{Z_star}{Change-in-score anomaly scores}
#' @export
compute_z_scores <- function(spectra, s = 5, l = 10) {
  T <- length(spectra)
  Z <- rep(0, T)

  for (t in (l + 1):T) {
    short_vec <- get_context_vector(spectra[(t - s):(t - 1)])
    long_vec  <- get_context_vector(spectra[(t - l):(t - 1)])
    cur_vec   <- spectra[[t]]
    Zs <- 1 - sum(cur_vec * short_vec)
    Zl <- 1 - sum(cur_vec * long_vec)
    Z[t] <- max(Zs, Zl)
  }

  Z_star <- rep(0, T)

  #Skip the first score
  Z_star[(l + 2):T] <- pmax(Z[(l + 2):T] - Z[(l + 1):(T - 1)], 0)

  return(list(Z = Z, Z_star = Z_star))
}


#' Plot LAD Z and Z* Scores Over Time
#'
#' @param Z Vector of raw anomaly scores (Z)
#' @param Z_star Vector of Z* (difference in scores)
#' @param changepoints Optional vector of changepoint indices to annotate
#' @param title Optional plot title
#'
#' @return No return value; generates plots
#' @export
plot_lad_scores <- function(Z, Z_star, changepoints = NULL, title = "LAD Scores", skip_t=NULL) {
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  if (!is.null(skip_t) && skip_t > 1) {
    Z_plot <- Z
    Z_plot[1:(skip_t - 1)] <- NA
  } else {
    Z_plot <- Z
  }

  plot(Z_plot, type = "l", col = "blue",
       main = paste(title, "- Z Score"),
       ylab = "Z", xlab = "Time",
       ylim = range(Z_plot, na.rm = TRUE))  # zooms to min/max of Z

  if (!is.null(changepoints)) abline(v = changepoints, col = "red", lty = 2)

  plot(Z_star, type = "l", col = "darkgreen", main = paste(title, "- Z* Score"),
       ylab = "Z*", xlab = "Time")
  if (!is.null(changepoints)) abline(v = changepoints, col = "red", lty = 2)
}

#' Run Full LAD Analysis and Plotting
#'
#' @param adjacency_list A list of adjacency matrices (one per time step)
#' @param s Short-term sliding window size (default: 5)
#' @param l Long-term sliding window size (default: 10)
#' @param changepoints Optional: vector of known change point indices for annotation
#' @param title Optional title for the plots
#' @param k Optional: number of singular values to retain (NULL = full spectrum)
#'
#' @return A list with Z and Z* scores
#' @export
run_lad_analysis <- function(adjacency_list,
                             s = 5,
                             l = 10,
                             changepoints = NULL,
                             title = "LAD Analysis",
                             k = NULL) {
  spectra <- lapply(adjacency_list, get_laplacian_spectrum, k = k)
  scores <- compute_z_scores(spectra, s = s, l = l)
  plot_lad_scores(scores$Z, scores$Z_star, changepoints = changepoints, title = title, skip_t=12)
  return(scores)
}

