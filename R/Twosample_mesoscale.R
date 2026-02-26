#' Mesoscale testing
#'
#' Asymptotic tests for network two-sample testing, described in \href{https://arxiv.org/abs/2410.17046}{MacDonald et al., (2024)}.
#' Used to test a mesoscale null hypothesis for a prespecified collection
#' (hypothesis set) of edges, under an exponential family edge model.
#' It is compatible with either weighted (Gaussian) edge or binary edge networks
#' (specified with argument \code{edge_type}), undirected or directed networks
#' (\code{directed}), and networks with or without self loops (\code{self_loops}).
#' The mesoscale projection testing methodology uses edge variables outside the hypothesis set
#' to learn a projection of the edges in the hypothesis set, to reduce the dimension of the
#' test and improve test power. This test supports one network per sample.
#'
#' @param A A matrix or list of matrices, adjacency matrix(es) for the first sample.
#' @param B A matrix or list of matrices, adjacency matrix(es) for the second sample.
#' @param hyp_set the hypothesis set of node pairs, specified as one of: (1) a 2-element list
#' containing the incident row and column indices (rectangle), (2) a list of 2-element lists each containing
#' incident row and column indices (collection of rectangles), or (3) a 2-column matrix of (row,column) index
#' pairs (unstructured hypothesis set).
#' @param edge_type a string, either \code{'weighted'} or \code{'binary'} depending on the possible
#' edge values. The test with weighted edges is performed assuming Gaussian edges. Defaults to \code{'weighted'}.
#' @param d a positive integer, specifies the dimension of projection in the mesoscale test. Defaults to \code{1}.
#' @param directed a Boolean, specifies whether the network is directed. If \code{directed=FALSE}, the
#' testing methodology will be adjusted to account for symmetry. Defaults to \code{TRUE}.
#' @param self_loops a Boolean, if \code{FALSE} the test will ignore diagonal entries. Defaults to \code{TRUE}.
#' @param proj_type a string, either \code{'impute'} or \code{'one_step'} to choose the projection learning
#' approach. Recommendation is \code{'one_step'} for rectangular hypothesis sets, \code{'impute'} for unstructured
#' hypothesis sets. Defaults to \code{'impute'}.
#' @param var_type a string, \code{'basic'}, \code{'orth'}, or \code{'quasi'} to choose the variance estimation approach. For
#' binary edges, \code{'quasi'} corresponds to fitting with a quasibinomial GLM.
#' For weighted edges, \code{'orth'} calibrates the variance based on the variability of the edges orthogonal to the projected data, while
#' \code{'basic'} uses the variability within the projected data. Defaults to \code{'basic'}
#' @param centered a Boolean, specifies whether to center the projections to ignore differences in the overall
#' edge mean on the mesoscale set. Defaults to \code{FALSE}.
#' @param masked_set a masked set of node pairs ignored in both the projection learning and testing phases
#' of the mesoscale methodology. Specified in the same format as \code{hyp_set}. Defaults to the empty set.
#'
#' @return A list containing:
#' \item{p.value}{The \eqn{p}-value for the test}
#' \item{statistic}{The asymptotically chi-square test statistic}
#' \item{df}{The degrees of freedom of the chi-square null distribution}
#'
#' @export
#'
#' @examples
#' model <- list(name='SBM',B=matrix(c(0.5,0.3,0.3,0.5),2,2),Pi=c(0.5,0.5),directed=TRUE,self_loops=TRUE)
#' data <- Simulate_netmodel(n=10,model,m=4,twosample=TRUE)
#' A <- data$A1; B <- data$A2
#'
#' # hypothesis set specified as a rectangle
#' meso1 <- Mesoscale_test(A,B,
#'                         hyp_set=list(1:4,1:4),
#'                         edge_type='binary',d=1,directed=TRUE,self_loops=TRUE,
#'                         proj_type='one_step',
#'                         var_type='basic')
#'
#' # same hypothesis set specified as unstructured indices
#' meso2 <- Mesoscale_test(A,B,
#'                         hyp_set=cbind(c(1,1,1,2,2,3),c(2,3,4,3,4,4)),
#'                         edge_type='binary',d=1,directed=TRUE,self_loops=TRUE,
#'                         proj_type='one_step',
#'                         var_type='basic')
#' # check that meso1 and meso2 give the same result
#'
#' # hypothesis set specified as two rectangles, test with centering
#' meso3 <- Mesoscale_test(A,B,
#'                         hyp_set=list(list(1:2,1:2),list(9:10,9:10)),
#'                         edge_type='binary',d=1,directed=TRUE,self_loops=FALSE,
#'                         proj_type='impute',
#'                         var_type='quasi',
#'                         centered=TRUE)
#'
Mesoscale_test <- function(A,B,
                           hyp_set, # a rectangle (list of {row,col}), a 2-column matrix of entries, or a list if rectangles
                           edge_type='weighted',# or binary
                           # other options
                           d=1,
                           directed = TRUE,
                           self_loops = TRUE, # overrides diagonal entries in the hyp_set, masks them in the imputation routine
                           proj_type='impute',# or onestep
                           var_type='basic', # or 'quasi'
                           centered=FALSE,
                           masked_set=list(NULL,NULL) # a rectangle (list of {row,col}), a 2-column matrix of entries, or a list if rectangles
){
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  B <- checklist(B)
  # convert hypothesis set and masked set of node pairs to two-column matrix format
  nodepair_indices <- hyp_set_to_indices(hyp_set,masked_set,directed,self_loops)
  hyp_indices <- nodepair_indices$hyp_indices
  masked_indices <- nodepair_indices$masked_indices
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # group means with transformation for binary edges
  if(edge_type=='binary'){
    Abar <- logit(apply(list_to_array(A),c(1,2),bmean))
    Bbar <- logit(apply(list_to_array(B),c(1,2),bmean))
  }
  else{
    Abar <- Reduce('+',A)/m1
    Bbar <- Reduce('+',B)/m2
  }
  # learn left and right-hand side projections
  if(proj_type=='impute'){
    hyp_proj <- Subspace_impute(Abar,Bbar,d,
                                hyp_indices,
                                masked_indices,
                                self_loops,
                                directed,
                                centered)
  }
  else{
    hyp_proj <- Subspace_onestep(Abar,Bbar,d,
                                 hyp_indices,
                                 masked_indices,
                                 self_loops,
                                 directed,
                                 centered)
  }
  # test procedures
  if(edge_type=='weighted'){
    test_out <- Weight_meso(A,B,
                            hyp_indices,
                            hyp_proj,
                            var_type,
                            directed,
                            centered)
  }
  else{
    test_out <- Binary_meso(A,B,
                            hyp_indices,
                            hyp_proj,
                            var_type,
                            directed,
                            centered)
  }
  # returns acceptance/rejection decision and a pvalue
  out <- list()
  out$p.value <- test_out$pval
  out$statistic <- test_out$stat
  out$df <- test_out$dfs
  return(out)
}
