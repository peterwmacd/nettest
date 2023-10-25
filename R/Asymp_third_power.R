Asymp_third_power <- function(A, B, sig, r = 2) {

  # size for each of the group
  m1 = length(A) # m1
  m2 = length(B) # m2


  # Compute sample average for each group
  Auij = Reduce("+", A) / length(A)
  Buij = Reduce("+", B) / length(B)

  # spectral clustering
  # of Auij and Buij based on the index



  # get n
  n = nrow(Auij) # nrow(Auij) should be = nrow(Buij)

  # Calculate Zij
  Zij = matrix(0, n, n)

  ## Calculate Zij, When i!=j, non-diag
  numer = Auij - Buij # {A_bar}1, ij - {A_bar}2, ij
  denom = sqrt(n * ((Auij) * (1 - Auij) / m1 + Buij * (1 - Buij) / m2))
  nondiag_z = numer / denom

  ## Calculate Zij, when i = j, diag
  Bij = matrix(0, n, n)
  for (i in 1:n) {
    if (runif(1) < 0.5) {
      Bij[i, i] = -1/sqrt(n)
    } else {
      Bij[i, i] = 1/sqrt(n)
    }
  }

  # Get Zij
  Zij = nondiag_z
  diag(Zij) = diag(Bij)
  print(Zij)




  # Compute test statistics
  Z_cubed <- Zij %*% Zij %*% Zij
  test.stat = sum(diag(Z_cubed)) / sqrt(15) # trace operator: sum of the diagonal elements
  # p value
  p.val <- 1 - pnorm(test.stat)
  test <- ifelse(p.val <= sig, 1, 0)

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

