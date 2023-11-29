

Asymp_third_power <- function(A, B, sig, r = 2, sig = 0.05, pmetric_method = "nbd") {

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
    if (runif(1) < 0.5) {
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
  p.val <- 1 - pnorm(test.stat)
  test <- ifelse(p.val <= sig, 1, 0)


  # # spectral clustering
  # # of Auij and Buij based on the index
  # idx = spectral_clus(list(a), list(b), 2)
  # Zij = matrix(0, 4, 4)

  # for (i in 1:2){
  #   for (j in 1:2) {
  #     if (i != j){
  #       curr_a = a[idx == i, idx == j]
  #       curr_b = b[idx == i, idx == j]
  #       # Pij = mean(curr)
  #       #
  #       # curr = B[[1]][idx == i, idx == j]
  #       # Qij = mean(curr)
  #       #在分母的计算上，用 mean（curra）
  #       numer = mean(curr_a) - mean(curr_b) # {A_bar}1, ij - {A_bar}2, ij
  #       denom = sqrt(4 * ((curr_a) * (1 - curr_a) / 3 + curr_b * (1 - curr_b) / 3))
  #       nondiag_z = numer / denom
  #       Zij[idx == i, idx == j] = nondiag_z
  #     }
  #
  #     else{
  #       # curr = A[[1]][idx == i, idx == j]
  #       # Pij = sum(curr) / (dim(curr)[1] * (dim(curr)[1] - 1)) # Pij might = x / 0
  #       #
  #       # curr = B[[1]][idx == i, idx == j]
  #       # Qij = sum(curr) / (dim(curr)[1] * (dim(curr)[1] - 1))
  #       Bij = matrix(0, 4, 4)
  #       for (k in 1:4) {
  #         if (runif(1) < 0.5) {
  #           Bij[k, k] = -1/sqrt(4)
  #         } else {
  #           Bij[k, k] = 1/sqrt(4)
  #         }
  #       }
  #       # diag(Zij[idx == i, idx == j]) = diag(Bij)
  #       curr_bij = Bij[idx == i, idx == j]
  #       Zij[idx == i, idx == j] = curr_bij
  #
  #     }
  #   }
  # }


  return(c(test, p.val))
}






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

