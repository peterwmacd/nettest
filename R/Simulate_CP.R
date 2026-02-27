#' Simulation for DCSBM Network Sequences with a Changepoint (or Anomaly)
#'
#' Simulates a lists of adjacency matrices on \code{n} nodes from the DCSBM with binary edges.
#' Sequences can be specified with a changepoint and two parameter sets will be constructed: model 1 (before change)
#' and model 2 (after change).
#' Additional model and changepoint information is specified as part of the \code{CP} argument, see below.
#'
#' @param n An integer, number of nodes in each network.
#' @param m An integer, number of networks in the sequence.
#' @param model A named list, provides additional model information:
#' \describe{
#'     \item{type}{Type of changepoint, one of \code{'none'}, \code{'B'}, \code{'theta'}, \code{'Z'}, \code{'merge'},
#'     \code{'split'}, or \code{'kidneyegg'}. \code{B}, \code{Z}, and \code{theta} denote the parameter which changes
#'     between model 1 and model 2. \code{merge} is a special case of \code{Z}, which merges the first two
#'     communities into one. \code{split} is a special case of \code{Z}, which randomly splits the first community into
#'     two. \code{kidneyegg} is a special case of \code{B} with 2 communities, which changes the connection rate within
#'     community 1.}
#'     \item{t_start}{Time of change from model 1 to model 2.}
#'     \item{t_end}{Time immediately preceding reversion from model 2 back to model 1. Defaults to \code{m}, i.e. no reversion.
#'     Note that \code{CP$t_end=CP$t_start} will produce an anomaly for just one network in the sequence.}
#'     \item{directed}{A Boolean, are the network edges directed? Defaults to \code{FALSE}.}
#'     \item{self_loops}{A Boolean, are self-loops allowed? Defaults to \code{FALSE}.}
#'     \item{Z}{For \code{none}, \code{B} or \code{theta}, community membership matrix or \eqn{n \times d} under both
#'     model 1 and model 2. If unspecified, community memberships are generated using \code{model$Pi}.}
#'     \item{Z1}{For \code{Z}, community membership matrix or \eqn{n \times d} under model 1.
#'     If unspecified, community memberships are generated using \code{model$Pi1}.}
#'     \item{Z2}{For \code{Z}, community membership matrix or \eqn{n \times d} under model 2.
#'     If unspecified, community memberships are generated using \code{model$Pi2}.}
#'     \item{Pi}{\eqn{K}-vector of community membership probabilities. For \code{merge} and \code{split}, used to generate
#'     community memberships under model 1.}
#'     \item{Pi1}{\eqn{K}-vector of community membership probabilities. For \code{Z}, used to generate
#'     community memberships under model 1 if \code{Z1} is unspecified.}
#'     \item{Pi2}{\eqn{K}-vector of community membership probabilities. For \code{Z}, used to generate
#'     community memberships under model 2 if \code{Z2} is unspecified.}
#'     \item{B}{For \code{none}, \code{theta} or \code{Z}, \eqn{K \times K} matrix of connection probabilties under
#'     both model 1 and model 2. For \code{merge}, matrix of connection probabilities under model 1; connections under model 2 are
#'     computed by averaging the connections for the first two communities.}
#'     \item{B1}{For \code{B} or \code{split}, \eqn{K \times K} matrix of connection probabilties under model 1.}
#'     \item{B2}{For \code{B} or \code{split}, \eqn{K \times K} matrix of connection probabilties under model 2.}
#'     \item{theta}{For all types besides \code{theta}, \eqn{n}-vector of degree correction parameters under both model 1 and model 2.
#'     Defaults to a vector of \code{1}'s.}
#'     \item{theta1}{For \code{theta}, \eqn{n}-vector of degree correction parameters under model 1.
#'     Defaults to a vector of \code{1}'s.}
#'     \item{theta2}{For \code{theta}, \eqn{n}-vector of degree correction parameters under model 2.
#'     Defaults to a vector of \code{1}'s.}
#'     \item{n_egg}{For \code{kidneyegg}, number of nodes to be assigned to community 1.}
#'     \item{p}{For \code{kidneyegg}, baseline connection probability for all edges besides within community 1 under model 2.}
#'     \item{q}{For \code{kidneyegg}, connection probability within community 1 under model 2.}
#' }
#'
#' @return A list containing:
#' \item{A}{A list of \code{m} matrices, simulated adjacency matrix(es).}
#' \item{Z1}{An \eqn{n \times K} community membership matrix for model 1.}
#' \item{Z2}{An \eqn{n \times K} community membership matrix for model 2.}
#'
#' @export
#'
#' @examples
#' set.seed(12)
#' n <- 100; m <- 10
#' B <- matrix(0.3,2,2) + diag(0.2,2,2)
#' Bnew <- B; Bnew[2,2] <- 0.9
#' Bext <- matrix(0.3,3,3) + diag(0.2,3,3)
#' Pi <- rep(0.5,2)
#' theta <- runif(n,0.8,1.2)
#' thetanew <- theta + 0.1
#'
#' # (1) none
#' CP1 <- list(type='none',t_start=5,B=B,Pi=Pi)
#' seq1 <- Simulate_netCP(n,m,CP=CP1)
#' #plot(sapply(seq1$A,mean))
#'
#' # (2) B
#' CP2 <- list(type='B',t_start=5,B1=B,B2=Bnew,Pi=Pi)
#' seq2 <- Simulate_netCP(n,m,CP=CP2)
#' #plot(sapply(seq2$A,mean))
#'
#' # (3) theta
#' CP3 <- list(type='theta',t_start=5,B=B,Pi=Pi,theta1=theta,theta2=thetanew)
#' seq3 <- Simulate_netCP(n,m,CP=CP3)
#' #plot(sapply(seq3$A,mean))
#'
#' # (4) Z
#' CP4 <- list(type='Z',t_start=5,B=B,Pi1=Pi,Pi2=Pi)
#' seq4 <- Simulate_netCP(n,m,CP=CP4)
#' #plot(sapply(seq4$A,function(x){mean(x[Z_to_C(seq4$Z1)==1,Z_to_C(seq4$Z1)==1])}))
#'
#' # (5) merge
#' CP5 <- list(type='merge',t_start=5,B=B,Pi=Pi)
#' seq5 <- Simulate_netCP(n,m,CP=CP5)
#' #plot(sapply(seq5$A,function(x){mean(x[Z_to_C(seq5$Z1)==1,Z_to_C(seq5$Z1)==1])}))
#'
#' # (6) split
#' CP6 <- list(type='split',t_start=5,B1=B,B2=Bext,Pi=Pi)
#' seq6 <- Simulate_netCP(n,m,CP=CP6)
#' #plot(sapply(seq6$A,function(x){mean(x[Z_to_C(seq6$Z2)==1,Z_to_C(seq6$Z2)==2])}))
#'
#' # (7) kidneyegg anomaly
#' CP7 <- list(type='kidneyegg',t_start=5,t_end=5,n_egg=20,p=0.5,q=0.8)
#' seq7 <- Simulate_netCP(n,m,CP=CP7)
#' #plot(sapply(seq7$A,function(x){mean(x[Z_to_C(seq7$Z1)==1,Z_to_C(seq7$Z1)==1])}))
Simulate_CP <- function(n,m,CP=list()){
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
           m1$B <- CP$B
           m2$B <- B_merge(CP$B) # average connection over communities 1,2
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
# B <- matrix(0.3,2,2) + diag(0.2,2,2)
# Bnew <- B; Bnew[2,2] <- 0.9
# Bext <- matrix(0.3,3,3) + diag(0.2,3,3)
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
