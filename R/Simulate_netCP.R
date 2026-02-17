# Simulating from binary DCSBM sequences (discrete time) with one changepoint or one anomaly, 
# wrapping around Simulate_netmodel
# all require n = #nodes, m = length of sequence, optional self_loops, directed 
# both default to FALSE (part of CP)
# t_start, t_end (default to m) for changepoint sequence, t_end=t_start gives a 
# sequence with a single anomaly

# CP types:
# - 'none': default, no CP/anomaly
# - 'B': specify B1,B2
# - 'theta': specify theta1,theta2
# - 'Z': specify Z1,Z2 
# - 'merge': special case of membership, specify B,Pi generate and merge communities 1+2
# - 'split': special case of membership, specify B,Pi, generate and randomly split community 1
# - 'kidneyegg': special case of 'B', provide (p,q,n_egg)

# main function
Simulate_netCP <- function(n,m,CP=list()){
  # default time_end = m
  if(is.null(CP$t_end)){
    CP$t_end <- m
  }
  # default undirected
  if(is.null(CP$directed)){
    CP$directed <- FALSE
  }
  # default no self loops
  if(is.null(CP$self_loops)){
    CP$self_loops <- FALSE
  }
  # initialize model info
  m1 <- m2 <- list(name='SBM',directed=CP$directed,self_loops=CP$self_loops)
  # CP switch
  switch(CP$type,
         none = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate memberships if Z is unspecified
           if(is.null(CP$Z)){
             K <- length(CP$Pi)
             C <- sample(1:K,n,replace=TRUE,prob=CP$Pi)
             m1$Z <- m2$Z <- C_to_Z(C,K)
           }
           else{
             m1$Z <- m2$Z <- CP$Z
           }
           # populate additional model info
           m1$B <- m2$B <- CP$B
         },
         B = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate memberships if Z is unspecified
           if(is.null(CP$Z)){
             K <- length(CP$Pi)
             C <- sample(1:K,n,replace=TRUE,prob=CP$Pi)
             m1$Z <- m2$Z <- C_to_Z(C,K)
           }
           else{
             m1$Z <- m2$Z <- CP$Z
           }
           # populate additional model info
           m1$B <- CP$B1; m2$B <- CP$B2
         },
         theta = {
           # generate memberships if Z is unspecified
           if(is.null(CP$Z)){
             K <- length(CP$Pi)
             C <- sample(1:K,n,replace=TRUE,prob=CP$Pi)
             m1$Z <- m2$Z <- C_to_Z(C,K)
           }
           else{
             m1$Z <- m2$Z <- CP$Z
           }
           # populate additional model info
           m1$B <- m2$B <- CP$B
           m1$theta <- CP$theta1; m2$theta <- CP$theta2
         },
         Z = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate memberships if Z1 is unspecified
           if(is.null(CP$Z1)){
             K1 <- length(CP$Pi1)
             C1 <- sample(1:K1,n,replace=TRUE,prob=CP$Pi1)
             m1$Z <- C_to_Z(C1,K1)
           }
           else{
             m1$Z <- CP$Z1
           }
           # generate memberships if Z2 is unspecified
           if(is.null(CP$Z2)){
             K2 <- length(CP$Pi2)
             C2 <- sample(1:K2,n,replace=TRUE,prob=CP$Pi2)
             m2$Z <- C_to_Z(C2,K2)
           }
           else{
             m2$Z <- CP$Z2
           }
           # populate additional model info
           m1$B <- m2$B <- CP$B
         },
         merge = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate Z1 from Pi and merge
           K <- length(CP$Pi)
           C <- sample(1:K,n,replace=TRUE,prob=CP$Pi)
           m1$Z <- C_to_Z(C,K)
           m2$Z <- Z_merge(m1$Z) # helper to merge communities 1,2
           # populate additional model info
           if(is.null(CP$B2)){
            m1$B <- CP$B
            m2$B <- B_merge(CP$B) # average connection over communities 1,2
           }
           else{
             m1$B <- CP$B1; m2$B <- CP$B2
           }
         },
         split = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate Z2 from Pi and merge
           Pi2 <- Pi_split(CP$Pi)
           K <- length(Pi2)
           C <- sample(1:K,n,replace=TRUE,prob=Pi2)
           m2$Z <- C_to_Z(C,K)
           m1$Z <- Z_merge(m2$Z) # helper to merge communities 1,2
           # populate additional model info
            m1$B <- CP$B1; m2$B <- CP$B2
         },
         kidneyegg = {
           # populate default degree parameters
           if(is.null(CP$theta)){
             m1$theta <- m2$theta <- rep(1,n)
           }
           else{
             m1$theta <- m2$theta <- CP$theta
           }
           # generate Z from n_egg
           Z <- matrix(0,n,2)
           Z[1:CP$n_egg,1] <- 1; Z[(CP$n_egg+1):n,2] <- 1
           m1$Z <- m2$Z <- Z
           # generate B1, B2 from p,q
           m1$B <- m2$B <- matrix(CP$p,2,2)
           m2$B[1,1] <- CP$q # chatter among egg-egg edges
         },
         {
           stop("CP type unavailable.")
         })
  # generate adjacency matrix sample(s)
  A1 <- checklist(Simulate_netmodel(n,m1,CP$t_start-1)$A)
  A2 <- checklist(Simulate_netmodel(n,m2,CP$t_end-CP$t_start+1)$A)
  if(CP$t_end == m){
    A3 <- NULL
  }
  else{
    A3 <- checklist(Simulate_netmodel(n,m1,m-CP$t_end)$A)
  }
  # return adjacency and community memberships (before/after CP)
  return(list(A=c(A1,A2,A3),Z1=m1$Z,Z2=m2$Z))
}

# #### test instances ####
# set.seed(12)
# # all with default undirected, no self loops
# # all with n=40 nodes
# n <- 100; m <- 10
# B <- B_balanced(2,0.5,0.3) # NOTE: uses helper B_balanced
# Bnew <- B; Bnew[2,2] <- 0.9
# Bext <- B_balanced(3,0.5,0.3)
# Pi <- rep(0.5,2)
# theta <- runif(n,0.8,1.2)
# thetanew <- theta + 0.1
# 
# # (1) none
# CP1 <- list(type='none',t_start=5,B=B,Pi=Pi)
# seq1 <- Simulate_netCP(n,m,CP=CP1)
# plot(sapply(seq1$A,mean))
# 
# # (2) B
# CP2 <- list(type='B',t_start=5,B1=B,B2=Bnew,Pi=Pi)
# seq2 <- Simulate_netCP(n,m,CP=CP2)
# plot(sapply(seq2$A,mean))
# 
# # (3) theta
# CP3 <- list(type='theta',t_start=5,t_end=6,B=B,Pi=Pi,theta1=theta,theta2=thetanew)
# seq3 <- Simulate_netCP(n,m,CP=CP3)
# plot(sapply(seq3$A,mean))
# 
# # (4) Z
# CP4 <- list(type='Z',t_start=5,B=B,Pi1=Pi,Pi2=Pi)
# seq4 <- Simulate_netCP(n,m,CP=CP4)
# plot(sapply(seq4$A,function(x){mean(x[Z_to_C(seq4$Z1)==1,Z_to_C(seq4$Z1)==1])}))
# 
# # (5) merge
# CP5 <- list(type='merge',t_start=5,B=B,Pi=Pi)
# seq5 <- Simulate_netCP(n,m,CP=CP5)
# plot(sapply(seq5$A,function(x){mean(x[Z_to_C(seq5$Z1)==1,Z_to_C(seq5$Z1)==1])}))
# 
# # (6) split
# CP6 <- list(type='split',t_start=5,B1=B,B2=Bext,Pi=Pi)
# seq6 <- Simulate_netCP(n,m,CP=CP6)
# plot(sapply(seq6$A,function(x){mean(x[Z_to_C(seq6$Z2)==1,Z_to_C(seq6$Z2)==2])}))
# 
# # (7) kidneyegg anomaly
# CP7 <- list(type='kidneyegg',t_start=5,t_end=5,n_egg=20,p=0.5,q=0.8)
# seq7 <- Simulate_netCP(n,m,CP=CP7)
# plot(sapply(seq7$A,function(x){mean(x[Z_to_C(seq7$Z1)==1,Z_to_C(seq7$Z1)==1])}))
