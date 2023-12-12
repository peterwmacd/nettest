# Asymptotic 3rd power test (Chen et al. 2020)
#' @export
Asymp_third_power <- function(A, B, sig=0.05, r = 2, pmetric_method = "nbd") {

  m1 = length(A) # m1
  m2 = length(B) # m2

  # spectral clustering using
  if (pmetric_method == "LG"){
    estr_1 = graphon::est.LG(A, r)
    estr_2 = graphon::est.LG(B, r) # still have NA

  } else if (pmetric_method == "SBA"){
    estr_1 = graphon::est.SBA(A, delta = 0.2)
    estr_2 = graphon::est.SBA(B, delta = 0.2)# all NA

  } else if (pmetric_method == "nbd"){
    estr_1 = graphon::est.nbdsmooth(A)
    estr_2 = graphon::est.nbdsmooth(B)# good, no NA

  } else if (pmetric_method == "USVT"){
    estr_1 = graphon::est.USVT(A, eta = 0.001)
    estr_2 = graphon::est.USVT(B, eta = 0.001)# all NA(0.001, 0.1, 0.9)

  } else if (pmetric_method == "spectral_clus"){
    # Added spectral clustering option w/existing helper function
    Abar <- Reduce('+',A) / length(A)
    Bbar <- Reduce('+',B) / length(B)
    # spectral clustering
    idx <- spectral_clus((Abar+Bbar)/2,r)
    # store estimates
    estr_1 = list(P=block_avg(Abar,idx))
    estr_2 = list(P=block_avg(Bbar,idx))
  } else {print("unexpected input for pmetric_method")}

  P_1 = estr_1$P
  P_2 = estr_2$P

  # Compute sample average for the two group
  A1ij = Reduce("+", A) / length(A)
  A2ij = Reduce("+", B) / length(B)

  Zij = matrix(0, dim(A1ij)[1], dim(A1ij)[1])

  num = A1ij - A2ij
  denom = sqrt(((P_1 %*% (1 - P_1)) / m1 + (P_2 %*% (1 - P_2)) / m2) * dim(A1ij)[1])

  Zij = num / denom

  Bij = matrix(0, dim(A1ij)[1], dim(A1ij)[1])
  for (k in 1:dim(A1ij)[1]) {
    if (stats::runif(1) < 0.5) {
      Bij[k, k] = -1/sqrt(dim(A1ij)[1])
    }
    else {
      Bij[k, k] = 1/sqrt(dim(A1ij)[1])
    }
  }

  diag(Zij) = diag(Bij)

  Z_cubed <- Zij %*% Zij %*% Zij
  test.stat = sum(diag(Z_cubed)) / sqrt(15) # trace operator: sum of the diagonal elements
  # p value
  p.val <- 1 - stats::pnorm(test.stat)
  test <- ifelse(p.val <= sig, 1, 0)
  return(c(test, p.val))
}

