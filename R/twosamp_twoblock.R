# Simulating stochastic block model networks with two blocks

# helper: random binary, undirected network with no self loops
rand_net <- function(P){
  n <- nrow(P)
  A <- matrix(0,n,n)
  A[lower.tri(P)] <- as.numeric(rbinom(n=sum(lower.tri(P)),
                                       size=1,
                                       prob=P[lower.tri(P)]))
  A[upper.tri(P)] <- t(A)[upper.tri(P)]
  return(A)
}

twosamp_twoblock <- function(n,m,p,q,epsilon){
  # block size
  n2 <- floor(n/2)
  # expected adjacency matrices
  # group 1
  P_G <- matrix(c(p,q,q,p),2,2) %x% matrix(1,n2,n2)
  # group 2
  P_H <- matrix(c(p,q,q,p+epsilon),2,2) %x% matrix(1,n2,n2)
  # generate networks
  # group 1
  A_G <- lapply(1:m,function(x){rand_net(P_G)})
  # group 2
  A_H <- lapply(1:m,function(x){rand_net(P_H)})
  return(list(A_G=A_G,A_H=A_H))
}
