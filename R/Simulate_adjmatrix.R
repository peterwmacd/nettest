# code to generate snapshots from different (independent edge) graph models
# after precalculating P = E(A)

# helper: simulate from general IER after constructing P = E(A)
Simulate_ier <- function(P,directed=FALSE,self_loops=FALSE){
  # dimension
  n <- nrow(P)
  if(directed){
    # sample A
    A <- matrix(as.numeric(stats::runif(n*n) < P),n,n)
    # remove self loops 
    if(!self_loops){
      diag(A) <- 0
    }
  }
  else{
    # extract upper triangle (w or w/o self loops)
    Ptri <- pracma::triu(P,k=ifelse(self_loops,0,1))
    # sample A (upper triangle)
    temp <- matrix(as.numeric(stats::runif(n*n) < Ptri),n,n)
    # symmetrize
    A <- temp + t(temp)
  }
  return(A)
}

# helper: special function for simulating a sample of overdispersed binary
# adjacency matrices: from general IER after constructing P = E(A)
Simulate_ier_od <- function(P,m,dispersion=1,
                            directed=FALSE,self_loops=FALSE){
  # dimension
  n <- nrow(P)
  # beta dispersion parameter
  betafac <- (m - dispersion)/(dispersion - 1)
  # 4d array to populate
  A <- array(0,c(n,n,m))
  for(ii in 1:n){
    if(!directed){
      for(jj in 1:ii){
        # generate beta success prob
        pr <- stats::rbeta(1,betafac*P[ii,jj],betafac*(1-P[ii,jj]))
        etot <- stats::rbinom(1,m,pr)
        eind <- sample(1:m,etot,replace=FALSE)
        A[ii,jj,][eind] <- 1
        # symmetric mirror
        A[jj,ii,][eind] <- 1
      }
    }
    else{
      for(jj in 1:n){
        # generate beta success prob
        pr <- stats::rbeta(1,betafac*P[ii,jj],betafac*(1-P[ii,jj]))
        etot <- stats::rbinom(1,m,pr)
        eind <- sample(1:m,etot,replace=FALSE)
        A[ii,jj,][eind] <- 1
      }
    }
  }
  # array to list and return
  Alist <- array_to_list(A,self_loops)
  return(Alist)
}

# helper: simulate from general IER after constructing P = E(A)
Simulate_ier_sample <- function(P,m=1,dispersion=1,
                                directed=FALSE,self_loops=FALSE){
  # toggle for multisample
  if(m==1){
    Simulate_ier(P,directed,self_loops)
  }
  else{
    # toggle for overdispersion (default=1)
    if(dispersion > 1){
      Simulate_ier_od(P,m,dispersion,directed,self_loops)
    }
    else{
      replicate(m,Simulate_ier(P,directed,self_loops),simplify=FALSE)
    }
  }
}

# helper: symmetric Gaussian noise with variance dispersion
sym_noise <- function(n,dispersion=1){
  nhalf <- n*(n+1)/2
  E <- matrix(0,n,n)
  E[upper.tri(E,diag=TRUE)] <- stats::rnorm(nhalf,sd=sqrt(dispersion))
  E[lower.tri(E,diag=FALSE)] <- t(E)[lower.tri(E,diag=FALSE)]
  return(E)
}

# helper: simulate from general Gaussian edge model after constructing E(A) = P,
# entry variance given by dispersion
Simulate_gaussnet <- function(P,dispersion=1,
                              directed=FALSE,self_loops=FALSE){
  # dimension
  n <- nrow(P)
  if(directed){
    # sample A
    A <- P + matrix(stats::rnorm(n^2,sd=sqrt(dispersion)),n,n)
    # remove self loops 
    if(!self_loops){
      diag(A) <- 0
    }
  }
  else{
    # extract upper triangle (w or w/o self loops)
    Ptri <- pracma::triu(P,k=ifelse(self_loops,0,1))
    # sample A (upper triangle)
    A <- P + sym_noise(n,dispersion)
    # remove self loops 
    if(!self_loops){
      diag(A) <- 0
    }
  }
  return(A)
}

# helper: simulate from general Gaussian edge model after constructing E(A) = P,
# entry variance given by dispersion
Simulate_gaussnet_sample <- function(P,m=1,dispersion=1,
                                     directed=FALSE,self_loops=FALSE){
  # toggle for multisample
  if(m==1){
    Simulate_gaussnet(P,dispersion,directed,self_loops)
  }
  else{
    replicate(m,Simulate_gaussnet(P,dispersion,directed,self_loops),simplify=FALSE)
  }
}