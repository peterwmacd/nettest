#' Get Laplacian Spectrum (with options for normalization and truncation)
#'
#' @description
#' Compute the eigenvalue spectrum of either the unnormalized Laplacian
#' (\eqn{L = D - A}) or the normalized Laplacian
#' (\eqn{L_{\text{sym}} = I - D^{-1/2} A D^{-1/2}}).
#' Supports truncation to the first/last \code{k} eigenvalues and optional
#' L2 normalization.
#'
#' @param A Numeric adjacency matrix (square, symmetric).
#' @param laplacian Character. Either \code{"normalized"} or \code{"unnormalized"}.
#' @param k Optional integer. Number of eigenvalues to keep (default = NULL = all).
#' @param which Character. Either \code{"smallest"} or \code{"largest"} when truncating.
#' @param normalize Logical. If TRUE, L2-normalize the spectrum vector (Huang default).
#'   If FALSE, return raw eigenvalues (Euclidean variant).
#'
#' @return Numeric vector of eigenvalues (possibly truncated, possibly normalized).
#' @export
get_laplacian_spectrum <- function(
    A,
    laplacian = c("normalized", "unnormalized"),
    k = NULL,
    which = c("smallest", "largest"),
    normalize = TRUE
) {
  laplacian <- match.arg(laplacian)
  which <- match.arg(which)
  A <- as.matrix(A)

  if (laplacian == "normalized") {
    deg <- rowSums(A)
    inv_sqrt_deg <- 1 / sqrt(pmax(deg, 1e-12))
    Dm12 <- diag(inv_sqrt_deg)
    L <- diag(nrow(A)) - Dm12 %*% A %*% Dm12
  } else {
    D <- diag(rowSums(A))
    L <- D - A
  }

  ev <- eigen(L, only.values = TRUE, symmetric = TRUE)$values
  ev <- sort(ev, decreasing = FALSE)

  if (!is.null(k)) {
    k <- min(k, length(ev))
    if (which == "smallest") {
      ev <- ev[seq_len(k)]
    } else {
      ev <- ev[(length(ev) - k + 1):length(ev)]
    }
  }

  if (isTRUE(normalize)) {
    nrm <- sqrt(sum(ev^2))
    if (nrm > 0) ev <- ev / nrm
  }

  return(ev)
}

#' Build a Context Vector from Spectra
#'
#' @description
#' Aggregate a set of spectral vectors into a single context vector,
#' either by mean or via truncated SVD (first left singular vector).
#'
#' @param spectra_list List of numeric vectors of equal length.
#' @param method Character. Either \code{"svd"} (first singular vector)
#'   or \code{"mean"} (componentwise average).
#'
#' @return Numeric vector of same length as each input spectrum.
#' @export
get_context_vector <- function(spectra_list, method = c("svd", "mean")) {
  method <- match.arg(method)
  mat <- do.call(cbind, spectra_list)
  if (method == "svd") {
    return(svd(mat)$u[, 1])
  } else {
    return(rowMeans(mat))
  }
}

#' Compute LAD Z and Z* Scores (Cosine or Euclidean)
#'
#' @description
#' Given a sequence of spectral signatures, compute anomaly scores:
#' \itemize{
#'   \item Z: raw anomaly score (distance from short-term or long-term context).
#'   \item Z*: change-in-score anomaly (positive increments in Z).
#' }
#' Supports cosine dissimilarity (Huang default) or Euclidean distance.
#'
#' @param spectra List of spectral vectors (equal length).
#' @param s Integer. Short-term window size (default = 5).
#' @param baseline_window Integer vector. Indices used for long-term baseline.
#' @param distance Character. Either \code{"cosine"} or \code{"euclidean"}.
#' @param context_method Character. Either \code{"svd"} or \code{"mean"} for context building.
#'
#' @return List with components \code{Z} and \code{Z_star}.
#' @export
compute_z_scores <- function(
    spectra,
    s = 5,
    baseline_window = 1:25,
    distance = c("cosine", "euclidean"),
    context_method = c("svd", "mean")
) {
  distance <- match.arg(distance)
  context_method <- match.arg(context_method)

  Tn <- length(spectra)
  Z <- rep(NA_real_, Tn)

  long_vec <- get_context_vector(spectra[baseline_window], method = context_method)

  cos_dis <- function(a, b) {
    denom <- sqrt(sum(a*a)) * sqrt(sum(b*b))
    if (denom == 0) return(0)
    1 - abs(sum(a * b) / denom)
  }
  euc_dis <- function(a, b) sqrt(sum((a - b)^2))

  start_t <- max(baseline_window) + s + 1

  if (start_t <= Tn) {
    for (t in start_t:Tn) {
      short_vec <- get_context_vector(spectra[(t - s):(t - 1)], method = context_method)
      cur_vec   <- spectra[[t]]

      Zs <- if (distance == "cosine") cos_dis(cur_vec, short_vec) else euc_dis(cur_vec, short_vec)
      Zl <- if (distance == "cosine") cos_dis(cur_vec, long_vec)   else euc_dis(cur_vec, long_vec)
      Z[t] <- max(Zs, Zl)
    }
  }

  Z_star <- rep(NA_real_, Tn)
  idx <- which(!is.na(Z))
  if (length(idx) >= 2) {
    Z_star[idx[-1]] <- pmax(Z[idx[-1]] - Z[idx[-length(idx)]], 0)
  }

  list(Z = Z, Z_star = Z_star)
}





#' Plot LAD Z and Z* with Shewhart limits
#'
#' @param Z numeric vector of Z scores
#' @param Z_star numeric vector of Z* scores
#' @param l long-window size (used to compute burn-in = max(baseline) + l)
#' @param baseline_window integer indices for baseline (e.g. 1:25)
#' @param changepoints optional integer vector to annotate with dashed v-lines
#' @param title base title stem (two panels will be titled " - Z" / " - Z*")
#' @export
plot_lad_scores <- function(Z, Z_star, l = 10, baseline_window = 1:15,
                            changepoints = NULL, title = "LAD") {
  burn_in <- max(baseline_window) + l


  # Baseline stats
  mean_Z  <- mean(Z[baseline_window], na.rm = TRUE)
  sd_Z    <- sd(Z[baseline_window],   na.rm = TRUE)
  UCL_Z   <- mean_Z + 3 * sd_Z
  LCL_Z   <- mean_Z - 3 * sd_Z

  mean_Zs <- mean(Z_star[baseline_window], na.rm = TRUE)
  sd_Zs   <- sd(Z_star[baseline_window],   na.rm = TRUE)
  UCL_Zs  <- mean_Zs + 3 * sd_Zs

  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  # ---- Z panel ----
  ylim_Z <- range(c(Z, UCL_Z, LCL_Z), na.rm = TRUE)
  plot(Z, type = "l", col = "blue", main = paste0(title, " - Z"),
       ylab = "Z", xlab = "Time", ylim = ylim_Z)

  # Shade burn-in region
  rect(xleft = 0, xright = burn_in, ybottom = par("usr")[3], ytop = par("usr")[4],
       col = rgb(0.9, 0.9, 0.9, 0.5), border = NA)

  lines(Z, col = "blue")
  abline(h = mean_Z, col = "blue")
  abline(h = c(UCL_Z, LCL_Z), col = "red", lty = 2)
  if (!is.null(changepoints)) abline(v = changepoints, col = "red", lty = 2)

  # ---- Z* panel ----
  ylim_Zs <- range(c(Z_star, UCL_Zs, 0), na.rm = TRUE)
  plot(Z_star, type = "l", col = "darkgreen",
       main = paste0(title, " - Z*"), ylab = "Z*", xlab = "Time", ylim = ylim_Zs)

  # Shade burn-in region
  rect(xleft = 0, xright = burn_in, ybottom = par("usr")[3], ytop = par("usr")[4],
       col = rgb(0.9, 0.9, 0.9, 0.5), border = NA)

  lines(Z_star, col = "darkgreen")
  abline(h = mean_Zs, col = "blue")
  abline(h = UCL_Zs, col = "red", lty = 2)
  if (!is.null(changepoints)) abline(v = changepoints, col = "red", lty = 2)

  par(mfrow = c(1, 1))
}




#' Run LAD Analysis (Flexible Variant)
#'
#' @description
#' Full LAD pipeline: extract Laplacian spectra, compute Z/Z* anomaly
#' scores, and plot control charts. Allows switching between
#' Huang’s default (normalized spectra + cosine distance) and
#' Euclidean variants (raw spectra, Euclidean distance, truncated eigenvalues).
#'
#' @param adjacency_list List of adjacency matrices (one per time step).
#' @param s Short-term window size (default = 5).
#' @param l Long-term burn-in padding (default = 10).
#' @param changepoints Optional integer vector of known changepoints for annotation.
#' @param title Character. Title prefix for plots.
#' @param laplacian Character. "normalized" or "unnormalized".
#' @param k Optional integer. Number of eigenvalues to keep (NULL = all).
#' @param which Character. "smallest" or "largest" eigenvalues to keep.
#' @param normalize Logical. If TRUE, L2-normalize (Huang). If FALSE, keep raw.
#' @param distance Character. "cosine" or "euclidean".
#' @param context_method Character. "svd" or "mean".
#' @param baseline_window Integer vector. Indices for baseline (default 1:25).
#'
#' @return List with components \code{Z} and \code{Z_star}.
#' @export
run_lad_analysis <- function(
    adjacency_list,
    s = 5,
    l = 10,
    changepoints = NULL,
    title = "LAD Analysis",
    laplacian = c("normalized", "unnormalized"),
    k = NULL,
    which = c("smallest", "largest"),
    normalize = TRUE,
    distance = c("cosine", "euclidean"),
    context_method = c("svd", "mean"),
    baseline_window = 1:25
) {
  laplacian <- match.arg(laplacian)
  which <- match.arg(which)
  distance <- match.arg(distance)
  context_method <- match.arg(context_method)

  spectra <- lapply(
    adjacency_list,
    get_laplacian_spectrum,
    laplacian = laplacian,
    k = k,
    which = which,
    normalize = normalize
  )

  scores <- compute_z_scores(
    spectra,
    s = s,
    baseline_window = baseline_window,
    distance = distance,
    context_method = context_method
  )

  plot_lad_scores(
    scores$Z, scores$Z_star,
    l = l,
    baseline_window = baseline_window,
    changepoints = changepoints,
    title = title
  )
  scores
}


