# mesoscale network generation

# from gaussian generation fn
# n, m, d, dispersion (overdisp for binary, sigma for gaussian, defaults to 1),
# row_indices, col_indices, signal, directed, gaussian/binary toggle, distance/IP toggle

# helper: convert array to list of aligned matrices
array_to_list <- function(A){
  m <- dim(A)[3]
  lapply(1:m,function(kk){A[,,kk]})
}

# helper convert list of aligned matrices to 3d array
list_to_array <- function(B){
  n <- nrow(B[[1]])
  m <- length(B)
  array(unlist(B),c(n,n,m))
}

# helper: hollowize a square matrix (set diagonal entries to zero)
hollowize <- function(A){
  A - diag(diag(A))
}

# helper: expit ### NOTE: already defined in utils.R
expit <- function(x){
  1/(1+exp(-x))
}

# helper: symmetric Gaussian noise with sd sigma
sym_noise <- function(n,sigma=1){
  nhalf <- n*(n+1)/2
  E <- matrix(0,n,n)
  E[upper.tri(E,diag=TRUE)] <- stats::rnorm(nhalf,sd=sigma)
  E[lower.tri(E,diag=FALSE)] <- t(E)[lower.tri(E,diag=FALSE)]
  return(E)
}

# helper: random binary, possibly undirected undirected network with no self loops
sym_binom <- function(P,m=1,sym=TRUE){
  n <- nrow(P)
  A <- matrix(0,n,n)
  if(sym){
    A[lower.tri(P,diag=TRUE)] <- as.numeric(stats::rbinom(n=sum(lower.tri(P,diag=TRUE)),
                                                   size=m,
                                                   prob=P[lower.tri(P,diag=TRUE)]))/m
    A[upper.tri(P,diag=FALSE)] <- t(A)[upper.tri(P,diag=FALSE)]
  }
  else{
    A[1:length(A)] <- as.numeric(rbinom(n=length(P),
                                        size=m,
                                        prob=c(P)))/m
  }
  return(A)
}

# generating mesoscale two-samples of networks (returns 2 lists of networks)
#' @export
genMesoscale <- function(n,m,d, # dimensions
                         row_indices=NULL,col_indices=NULL, # indices for mesoscale test set
                         model='gaussian', # 'gaussian'/'logistic' toggle
                         directed=TRUE, # directed/undirected edge toggle
                         self_loops=TRUE,
                         signal=0, # signal control (on difference in posns)
                         dispersion=1, # dispersion (edge sd or overdispersion param)
                         simfunc='IP'){ # similarity function for posns (inner prod or 'distance')
  # 1. Generate positions
  if(!directed){
    # LHS positions
    X1 <- matrix(stats::rnorm(n*d),n,d)
    X2 <- matrix(stats::rnorm(n*d),n,d)
    X2[union(row_indices,col_indices),] <- X1[union(row_indices,col_indices),] + matrix(stats::rnorm(length(union(row_indices,col_indices))*d,sd=signal/sqrt(d)),length(union(row_indices,col_indices)),d)
    # RHS positions (symmetric)
    Y1 <- X1
    Y2 <- X2
  }
  else{
    # LHS positions
    X1 <- matrix(stats::rnorm(n*d),n,d)
    X2 <- matrix(stats::rnorm(n*d),n,d)
    X2[row_indices,] <- X1[row_indices,] + matrix(stats::rnorm(length(row_indices)*d,sd=signal/sqrt(d)),length(row_indices),d)
    # RHS positions
    Y1 <- matrix(stats::rnorm(n*d),n,d)
    Y2 <- matrix(stats::rnorm(n*d),n,d)
    Y2[col_indices,] <- Y1[col_indices,] + matrix(stats::rnorm(length(col_indices)*d,sd=signal/sqrt(d)),length(col_indices),d)
  }

  # 2. Compute Theta matrix
  if(model=='logistic'){
    if(simfunc=='distance'){
      # edge expectations
      if(!directed){
        Theta1 <- expit(as.matrix(stats::dist(X1)))
        Theta2 <- expit(as.matrix(stats::dist(X2)))
      }
      else{
        Theta1 <- expit(as.matrix(stats::dist(rbind(X1,Y1)))[1:n,-(1:n)])
        Theta2 <- expit(as.matrix(stats::dist(rbind(X2,Y2)))[1:n,-(1:n)])
      }
    }
    else{
      # edge expectations
      Theta1 <- expit(tcrossprod(X1,Y1))
      Theta2 <- expit(tcrossprod(X2,Y2))
    }
  }
  else{
    if(simfunc=='distance'){
      # edge expectations
      if(!directed){
        Theta1 <- as.matrix(stats::dist(X1))
        Theta2 <- as.matrix(stats::dist(X2))
      }
      else{
        Theta1 <- as.matrix(stats::dist(rbind(X1,Y1)))[1:n,-(1:n)]
        Theta2 <- as.matrix(stats::dist(rbind(X2,Y2)))[1:n,-(1:n)]
      }
    }
    else{
      # edge expectations
      Theta1 <- tcrossprod(X1,Y1)
      Theta2 <- tcrossprod(X2,Y2)
    }
  }
  # update if self loops are not allowed, NOTE: hollow Theta1 -> no self loops under the logistic model;
  # need additional hollowization for gaussian model
  if(!self_loops){
    Theta1 <- hollowize(Theta1)
    Theta2 <- hollowize(Theta2)
  }

  # 3. Compute observed signal
  if(!is.null(row_indices) & !is.null(row_indices)){
    signal_obs <- sqrt(sum((Theta1[row_indices,col_indices] - Theta2[row_indices,col_indices])^2))
  }
  else{
    signal_obs <- NA
  }

  # 4. Generate networks
  if(model=='logistic'){
    if(dispersion>1){ # overdispersed case
      # populate networks
      # beta dispersion
      betafac <- (m - dispersion)/(dispersion - 1)
      # store Thetas together in an array
      Theta <- array(c(Theta1,Theta2),c(n,n,2))
      # 4d array to populate
      A <- array(0,c(n,n,m,2))
      for(ii in 1:n){
        if(!directed){
          for(jj in 1:ii){
            for(gg in 1:2){
              # generate beta success prob
              pr <- stats::rbeta(1,betafac*Theta[ii,jj,gg],betafac*(1-Theta[ii,jj,gg]))
              etot <- stats::rbinom(1,m,pr)
              eind <- sample(1:m,etot,replace=FALSE)
              A[ii,jj,,gg][eind] <- 1
              # symmetric mirror
              A[jj,ii,,gg][eind] <- 1
            }
          }
        }
        else{
          for(jj in 1:n){
            for(gg in 1:2){
              # generate beta success prob
              pr <- stats::rbeta(1,betafac*Theta[ii,jj,gg],betafac*(1-Theta[ii,jj,gg]))
              etot <- stats::rbinom(1,m,pr)
              eind <- sample(1:m,etot,replace=FALSE)
              A[ii,jj,,gg][eind] <- 1
            }
          }
        }
      }
      # split array to return each group
      A1 <- A[,,,1]
      A2 <- A[,,,2]
    }
    else{
      # populate networks
      A1 <- A2 <- array(NA,c(n,n,m))
      for(kk in 1:m){
        A1[,,kk] <- sym_binom(Theta1,sym=!directed)
        A2[,,kk] <- sym_binom(Theta2,sym=!directed)
      }
    }
  }
  else{
    # populate networks
    A1 <- A2 <- array(NA,c(n,n,m))
    for(kk in 1:m){
      if(!directed){
        A1[,,kk] <- Theta1 + sym_noise(n,dispersion)
        A2[,,kk] <- Theta2 + sym_noise(n,dispersion)
      }
      else{
        A1[,,kk] <- Theta1 + matrix(stats::rnorm(n^2,sd=dispersion),n,n)
        A2[,,kk] <- Theta2 + matrix(stats::rnorm(n^2,sd=dispersion),n,n)
      }
      if(!self_loops){
        A1[,,kk] <- hollowize(A1[,,kk])
        A2[,,kk] <- hollowize(A2[,,kk])
      }
    }
  }

  # 5. Return output
  return(list(A1=array_to_list(A1),A2=array_to_list(A2),
              Theta1=Theta1,Theta2=Theta2,
              X1=X1,X2=X2,
              Y1=Y1,Y2=Y2,
              signal_obs=signal_obs,
              row_indices=row_indices,col_indices=col_indices))
}

# # testing with simulation cases from paper
# set.seed(12)
#
# # null Gaussian IP
# test1 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='gaussian',
#                       directed=TRUE,
#                       signal=0,
#                       dispersion=sqrt(50),
#                       simfunc='IP')
#
# # non-null Gaussian IP
# test2 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='gaussian',
#                       directed=TRUE,
#                       signal=1/sqrt(2),
#                       dispersion=sqrt(50),
#                       simfunc='IP')
#
# # non-null Gaussian distance, undirected
# test3 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='gaussian',
#                       directed=FALSE,
#                       self_loops=FALSE,
#                       signal=2/sqrt(2),
#                       dispersion=sqrt(50),
#                       simfunc='distance')
#
# # null logistic IP
# test4 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='logistic',
#                       directed=TRUE,
#                       signal=0,
#                       dispersion=1,
#                       simfunc='IP')
#
# # non-null logistic IP
# test5 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='logistic',
#                       directed=TRUE,
#                       signal=0.5/sqrt(2),
#                       dispersion=1,
#                       simfunc='IP')
#
# # non-null logistic IP, overdispersed
# test6 <- genMesoscale(n=100,m=10,d=3,
#                       row_indices=1:20,col_indices=71:100,
#                       model='logistic',
#                       directed=TRUE,
#                       self_loops=FALSE,
#                       signal=0.5/sqrt(2),
#                       dispersion=2,
#                       simfunc='IP')

