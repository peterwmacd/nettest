# Simulating from different network models, wrapping around Simulate_adjmatrix
# all require n = #nodes, optional self_loops, directed both default to FALSE
# optional m for repeated sampling: if m > 1 returns a list of adjacency matrices
# from repeated samples; if twosample=TRUE splits into two list of size m

# ER
# requires p

# (DC)SBM
# requires (B,Z) or (B,Pi), optional theta vector

# LSM
# requires Z or d (default rnorm or rdirichlet positions), link (identity or logit,
# only for binary), similarity (ip or dist + optional intercept alpha for binary dist)

# Gaussian dispersion = sigma^2, binary dispersion can induce overdispersion in 
# samples of sive m > 1

# Function for binary edge networks
Simulate_netmodel <- function(n,model=list(), #snapshot model
                              m=1,twosample=FALSE){ #repeated sampling
  # set default parameters
  # default undirected
  if(is.null(model$directed)){
    model$directed <- FALSE
  }
  # default no self loops
  if(is.null(model$self_loops)){
    model$self_loops <- FALSE
  }
  # default dispersion/sigma^2 = 1
  if(is.null(model$dispersion)){
    model$dispersion <- 1
  }
  # determine binary or gaussian from model$name
  if(model$name %in% c('ER','SBM','LSM')){
    binary <- TRUE
  }
  else{
    binary <- FALSE
  }
  # model switch
  switch(model$name,
         ER = {
           # dummy Z for return statement
           Z <- NULL
           # compute full probability matrix
           P <- matrix(model$p,n,n)
         },
         SBM = {
           # default degree correction parameters ==1
           if(is.null(model$theta)){
             model$theta <- rep(1,n)
           }
           # generate memberships if Z is unspecified
           if(is.null(model$Z)){
             K <- length(model$Pi)
             C <- sample(1:K,n,replace=TRUE,prob=model$Pi)
             model$Z <- C_to_Z(C,K)
           }
           # compute full probability matrix (with prob. clipping)
           P <- pclip((model$Z %*% (tcrossprod(model$B, model$Z))) * (model$theta %o% model$theta),0)
         },
         LSM = {
           # generate memberships if Z is unspecified
           if(is.null(model$Z)){
             if(model$link == 'identity'){
               model$Z <- MCMCpack::rdirichlet(n,rep(1,model$d))
             }
             else{
               model$Z <- matrix(rnorm(n*model$d),nrow=n)
             }
           }
           # compute similarity matrix
           if(model$similarity=='dist'){
             if(is.null(model$alpha)){
               model$alpha <- 0
             }
             S <- as.matrix(model$alpha - dist(model$Z))
           }
           else{
             S <- tcrossprod(model$Z)
           }
           # compute prob. matrix through link function
           if(model$link == 'identity'){
             P <- pclip(S,0)
           }
           else{
             P <- expit(S)
           }
         },
         wSBM = {
           # default degree correction parameters ==1
           if(is.null(model$theta)){
             model$theta <- rep(1,n)
           }
           # generate memberships if Z is unspecified
           if(is.null(model$Z)){
             K <- length(model$Pi)
             C <- sample(1:K,n,replace=TRUE,prob=model$Pi)
             model$Z <- C_to_Z(C,K)
           }
           # compute expected adjacency matrix
           P <- model$Z %*% (tcrossprod(model$B, model$Z)) * (model$theta %o% model$theta)
         },
         wLSM = {
           # generate memberships if Z is unspecified
           if(is.null(model$Z)){
             Z <- matrix(rnorm(n*model$d),nrow=n)
           }
           # compute expected adjacency matrix (no link)
           if(model$similarity=='dist'){
             P <- as.matrix(dist(Z))
           }
           else{
             P <- tcrossprod(Z)
           }
         },
         {
           stop("Model name unavailable.")
         })
  # generate adjacency matrix sample(s)
  if(binary){
    if(twosample){
      A1 <- Simulate_ier_sample(P,m,model$dispersion,
                                directed=model$directed,self_loops=model$self_loops)
      A2 <- Simulate_ier_sample(P,m,model$dispersion,
                                directed=model$directed,self_loops=model$self_loops)
    }
    else{
      A <- Simulate_ier_sample(P,m,model$dispersion,
                               directed=model$directed,self_loops=model$self_loops)
    }
  }
  else{
    if(twosample){
      A1 <- Simulate_gaussnet_sample(P,m,model$dispersion,
                                     directed=model$directed,self_loops=model$self_loops)
      A2 <- Simulate_gaussnet_sample(P,m,model$dispersion,
                                     directed=model$directed,self_loops=model$self_loops)
    }
    else{
      A <- Simulate_gaussnet_sample(P,m,model$dispersion,
                                    directed=model$directed,self_loops=model$self_loops)
    }
  }
  # return adjacency and communities/positions
  if(twosample){
    return(list(A1=A1,A2=A2,Z=model$Z))
  }
  else{
    return(list(A=A,Z=model$Z))
  }
}

#### test instances ####
# all with default undirected, no self loops
# all with n=40 nodes
n <- 40

# (1) ER
m1 <- list(name='ER',p=0.3)
dat1 <- Simulate_netmodel(n,model=m1)

# (2) SBM fixed Z
Z2 <- C_to_Z(sample(1:2,n,replace=TRUE),2)
B2 <- matrix(c(0.7,0.3,0.3,0.7),2,2)
m2 <- list(name='SBM',Z=Z2,B=B2)
dat2 <- Simulate_netmodel_binary(n,model=m2)

# (3) DCSBM random Z
Pi3 <- c(0.5,0.5)
theta3 <- runif(n,min=0.8,1.2)
m3 <- list(name='SBM',B=B2,Pi=Pi3,theta=theta3)
dat3 <- Simulate_netmodel(n,model=m3)

# (4) binary RDPG fixed Z
Z4 <- MCMCpack::rdirichlet(n,rep(1,2))
m4 <- list(name='LSM',Z=Z4,link='identity',similarity='ip')
dat4 <- Simulate_netmodel(n,model=m4)

# (5) binary LSM random Z (dist sim + logit link)
m5 <- list(name='LSM',d=2,link='logit',similarity='dist',alpha=1)
dat5 <- Simulate_netmodel(n,model=m5)

# (6) Gaussian, block constant/SBM
m6 <- list(name='wSBM',B=3*B2,Pi=Pi3)
dat6 <- Simulate_netmodel(n,model=m6)

# (7) Gaussian, LSM random Z (ip sim)
m7 <- list(name='wLSM',d=3,dispersion=3,similarity='ip')
dat7 <- Simulate_netmodel(n,m7)


