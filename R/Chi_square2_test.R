library(Matrix)
library(irlba)


# Chi-square2 Test ---------------------------------------------------------
# Note: 
#   All graphs are assumed to be unweighted, undirected, and defined on a common vertex set. 
#   Potential commands: pchisq(..., n * (n - 1) / 2, lower.tail = FALSE)	

# Input:
#   A: cell array containing networks in 1st population; each cell is a sparse adjacency matrix
#   B: cell array containing networks in 2nd population, each cell being a sparse adjacency matrix
#   sig: significance level for acceptance of null hypothesis

# Output:
#   acceptance/rejection decision
#   p-value for the chi2 test

# Example: Asymp_chi2(simu$A_G, simu$A_H, 0.05) 

Asymp_chi2 <- function(A, B, sig) {
  
  m.a = length(A)
  m.b = length(B)
  n <- dim(A[[1]])[1] # nodes
  
  vec.A = array(0, dim = c(m.a, n * n))
  vec.B = array(0, dim = c(m.b, n * n))
  
  for (i in 1:m.a){
    curr = A[[i]]
    curr[lower.tri(curr, diag = TRUE)] = 0
    vec.A[i, ] = t(as.vector(curr))
  }
  for (i in 1:m.b){
    curr = B[[i]]
    curr[lower.tri(curr, diag = TRUE)] = 0
    vec.B[i, ] = t(as.vector(curr))
  }
  
  idx_non_zero = (colSums(vec.A + vec.B) != 0)
  vec.A = vec.A[, idx_non_zero]
  vec.B = vec.B[, idx_non_zero]
  
  mean.diff = colMeans(vec.A) - colMeans(vec.B)
  numer = mean.diff ^ 2
  denom = apply(vec.A, MARGIN = 2, var) / m.a + 
    apply(vec.B, MARGIN = 2, var) / m.b
  
  #denom = 0
  ind_nonzero_dem = (denom != 0)
  numer = numer[ind_nonzero_dem]
  denom = denom[ind_nonzero_dem]
  
  test.stat = sum(numer / denom)
  p.val = pchisq(test.stat, df = n * (n - 1) / 2, lower.tail = FALSE)
  test <- ifelse(p.val <= sig, 1, 0)
  return(c(test, p.val))
}