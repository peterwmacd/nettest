#' Simulation for Two-Sample Networks with Mesoscale Differences
#'
#' Simulates two lists of adjacency matrices on \code{n} nodes from latent space models with
#' binary edges (with logit link) or weighted (Gaussian) edges. Expected networks can differ on a specified (rectangular) set of edges.
#' Additional model information is specified as part of the \code{model} argument, see below.
#'
#' @param n An integer, number of nodes in each network.
#' @param model A named list, provides additional model information:
#' \describe{
#'     \item{name}{Either \code{'LSM'} or \code{'wLSM'}, describing the edge type.}
#'     \item{hyp_set}{A length two list containing a vector of row indices and a vector of column indices, used to
#'     construct the rectangular set of differential edges.}
#'     \item{directed}{A Boolean, are the network edges directed? Defaults to \code{FALSE}.}
#'     \item{self_loops}{A Boolean, are self-loops allowed? Defaults to \code{FALSE}.}
#'     \item{dispersion}{Edge dispersion parameter. For weighted networks, corresponds to edge variance.}
#'     \item{d}{For \code{LSM} and \code{wLSM}, integer latent space dimension.}
#'     \item{signal}{Perturbation standard deviation for latent positions. Each entry of a latent position
#'     incident to the hypothesis set rows or columns is perturbed by Gaussian noise with variance \code{signal^2/d}.}
#'     \item{similarity}{Either \code{'ip'} or \code{'dist'}, corresponding to inner product and Euclidean distance similarity functions.}
#'     \item{alpha}{For \code{LSM} with \code{model$similarity='dist'}, intercept term for edge probabilities.}
#' }
#' @param m An integer, number of networks in each sample.
#'
#' @return A list containing:
#' \item{A1}{A list of \code{m} matrices, simulated adjacency matrix(es) for sample 1.}
#' \item{A2}{A list of \code{m} matrices, simulated adjacency matrix(es) for sample 2.}
#' \item{Z1}{An \eqn{n \times d} latent position matrix used to generate sample 1.}
#' \item{Z2}{An \eqn{n \times d} latent position matrix used to generate sample 2.}
#' \item{Y1}{An \eqn{n \times d} right-hand side latent position matrix used to generate sample 1 (only if \code{model$directed=TRUE}).}
#' \item{Y2}{An \eqn{n \times d} right-hand side latent position matrix used to generate sample 2 (only if \code{model$directed=TRUE}).}
#'
#' @export
#'
#' @examples
#' set.seed(12)
#' n <- 100; m <- 10
#' H <- list(1:20,71:100)
#'
#' # (1) null Gaussian IP
#' model1 <- list(name='wLSM',hyp_set=H,d=3,signal=0,dispersion=1,similarity='ip')
#' data1 <- Simulate_mesoscale(n,model1,m)
#' #image(data1$A1[[1]] - data1$A2[[1]])
#'
#' # (2) non-null Gaussian IP
#' model2 <- list(name='wLSM',hyp_set=H,d=3,signal=1/sqrt(2),dispersion=1,similarity='ip')
#' data2 <- Simulate_mesoscale(n,model2,m)
#' #image(data2$A1[[1]] - data2$A2[[1]])
#'
#' # (3) null Gaussian distance, undirected
#' model3 <- list(name='wLSM',hyp_set=H,d=3,signal=0,dispersion=1,similarity='dist',directed=FALSE)
#' data3 <- Simulate_mesoscale(n,model3,m)
#' #image(data3$Z1 - data3$Z2)
#'
#' # (4) null logistic IP
#' model4 <- list(name='LSM',hyp_set=H,d=3,signal=0,similarity='ip')
#' data4 <- Simulate_mesoscale(n,model4,m)
#' #table(data4$A1[[1]] - data4$A2[[1]])
#'
#' # (5) non-null logistic IP
#' model5 <- list(name='LSM',hyp_set=H,d=3,signal=0.5/sqrt(2),similarity='ip')
#' data5 <- Simulate_mesoscale(n,model5,m)
#' #table(data5$A1[[1]] - data5$A2[[1]])
#'
#' # (6) non-null logistic IP, overdispersed
#' model6 <- list(name='LSM',hyp_set=H,d=3,signal=0.5/sqrt(2),dispersion=2,similarity='ip')
#' data6 <- Simulate_mesoscale(n,model6,m)
#' #table(data6$A1[[1]] - data6$A2[[1]])
Simulate_mesoscale <- function(n,model=list(),m){
  # check hypothesis set
  if(!(length(model$hyp_set)==2)){
    stop('Only supports rectangular hypothesis sets')
  }
  # set default parameters
  # default undirected
  if(is.null(model$directed)){
    model$directed <- TRUE
  }
  # default no self loops
  if(is.null(model$self_loops)){
    model$self_loops <- TRUE
  }
  # default dispersion/sigma^2 = 1
  if(is.null(model$dispersion)){
    model$dispersion <- 1
  }
  # default signal = 0
  if(is.null(model$signal)){
    model$signal <- 0
  }
  # default alpha (dist model intercept) = 0
  if(is.null(model$alpha)){
    model$alpha <- 0
  }
  # default link function (logit for LSM)
  if(model$name=='LSM'){
    model$link <- 'logit'
  }
  # copy models for two-samples
  model1 <- model2 <- model
  # E(A) for each sample
  # recover row/col indices, perturbation std.dev (psd)
  row_indices <- model$hyp_set[[1]]; col_indices <- model$hyp_set[[2]]
  all_indices <- union(row_indices,col_indices)
  psd <- model$signal/sqrt(model$d)
  # generate positions
  if(model$directed){
    # LHS positions
    Z1 <- matrix(stats::rnorm(n*d),n,d); Z2 <- matrix(stats::rnorm(n*d),n,d)
    Z2[row_indices,] <- perturb_mat(Z1[row_indices,],sd=psd)
    # RHS positions
    Y1 <- matrix(stats::rnorm(n*d),n,d); Y2 <- matrix(stats::rnorm(n*d),n,d)
    Y2[col_indices,] <- perturb_mat(Y1[col_indices,],sd=psd)
    # store positions
    model1$Z <- Z1; model2$Z <- Z2
    model1$Y <- Y1; model2$Y <- Y2
  }
  else{
    # LHS positions
    Z1 <- matrix(stats::rnorm(n*d),n,d); Z2 <- matrix(stats::rnorm(n*d),n,d)
    Z2[all_indices,] <- perturb_mat(Z1[all_indices,],sd=psd)
    # store positions
    model1$Z <- Z1; model2$Z <- Z2
  }
  # generate networks
  A1 <- Simulate_netmodel(n,model1,m)$A; A2 <- Simulate_netmodel(n,model2,m)$A
  # store and return
  if(model$directed){
  return(list(A1=A1,A2=A2,
              Z1=Z1,Z2=Z2,
              Y1=Y1,Y2=Y2))
  }
  else{
    return(list(A1=A1,A2=A2,
                Z1=Z1,Z2=Z2))
  }
}
