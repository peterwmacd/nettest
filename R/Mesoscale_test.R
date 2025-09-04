# mesoscale network testing

# binary_meso helper for logistic model
# takes the data, hyp_set, projections, variance option
# returns the pvalue and rejection decision
Binary_meso <- function(A,B,
                        hyp_indices,
                        hyp_proj,
                        var_type,
                        directed,
                        centered){
  # left/right objects
  Lproj <- hyp_proj$Lproj
  Rproj <- hyp_proj$Rproj
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # reset tilde if m is 1
  if(min(m1,m2)==1){
    var_type <- 'basic'
  }
  # hypothesis set entries as a matrix
  s_ind <- matrix(FALSE,n,n)
  s_ind[hyp_indices] <- TRUE
  # 'basic' design
  if(!directed){
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind & upper.tri(s_ind,diag=TRUE)
    # linear mapping and projection basis
    Gd <- Sym_span(s_ind,s_ind_data)
    UV <- pracma::orth(crossprod(Gd,(Rproj %x% Lproj)[c(s_ind),]))
    # need orth to remove possible rank deficiency
  }
  else{
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind
    # projection basis
    UV <- (Rproj %x% Lproj)[c(s_ind),]
  }
  # centering
  if(centered){
    UV <- cbind(rep(1/sqrt(nrow(UV)),nrow(UV)),center_orth(UV))
  }
  # dimensions
  ss <- nrow(UV)
  dd <- ncol(UV)
  # common subspace
  X <- rbind(UV,UV)
  # difference subspace
  W <- rbind(UV,-UV)
  # final design
  Xall <- cbind(X,W)

  # populate successes
  s_vecA <- rowSums(sapply(1:m1,function(kk){c(A[[kk]][s_ind_data],rep(0,ss))}))
  s_vecB <- rowSums(sapply(1:m2,function(kk){c(rep(0,ss),B[[kk]][s_ind_data])}))
  s_vec <- s_vecA + s_vecB
  # populate failures
  f_vec <- c(rep(m1,ss),rep(m2,ss)) - s_vec

  # fit logistic regression
  if(var_type=='quasi'){
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::quasibinomial())
  }
  else{
    moda <- stats::glm(cbind(s_vec,f_vec)~Xall-1,family=stats::binomial())
  }
  # get logistic regression coefficients
  gamma1 <- moda$coefficients[1:dd]
  if(centered){
    gamma2 <- moda$coefficients[-(1:(dd+1))]
  }
  else{
    gamma2 <- moda$coefficients[-(1:dd)]
  }
  # estimate covariance
  # linear model predictor
  lmhat <- c(UV %*% gamma1)
  # saturated model predictor
  prop_vec <- (s_vec[1:ss] + s_vec[-(1:ss)])/(m1+m2)
  prop_vec_reg <- sign(prop_vec - 0.5)*pmax(abs(prop_vec - 0.5) - (1/(2*(m1+m2))),0) + 0.5
  smhat <- Logit(prop_vec_reg)
  # covariances
  if(centered){
    Ghat <- 2*(t(UV[,-1,drop=FALSE]) %*% (diag(vexpit(lmhat)) %*% UV[,-1,drop=FALSE]))
    Fhat <- 2*(t(UV[,-1,drop=FALSE]) %*% (diag(vexpit(smhat)) %*% UV[,-1,drop=FALSE]))
  }
  else{
    Ghat <- 2*(t(UV) %*% (diag(vexpit(lmhat)) %*% UV))
    Fhat <- 2*(t(UV) %*% (diag(vexpit(smhat)) %*% UV))
  }
  # full inverse covariance
  iVhat <- Ghat %*% (solve(Fhat) %*% Ghat)
  # stat numerator
  num <- ((m1+m2)/2)*sum(gamma2 * (iVhat %*% gamma2))
  # overdispersion estimate as denominator
  den <- summary(moda)$dispersion
  # compile results and return
  out <- list()
  out$stat <- num/den
  if(centered){
    out$dfs <- dd-1
    out$pval <- stats::pchisq(out$stat,dd-1,lower.tail=FALSE)
  }
  else{
    out$dfs <- dd
    out$pval <- stats::pchisq(out$stat,dd,lower.tail=FALSE)
  }
  return(out)
}


# weight_meso2 helper for Gaussian model
# takes the data, hyp_set, projections, variance option (basic, quasi)
# returns the pvalue and rejection decision
## avoids using rectangular properties
## assumes n x d input projection matrices
## assumes hypothesis set passed as a 2 column matrix of indices
Weight_meso <- function(A,B,
                        hyp_indices,
                        hyp_proj,
                        var_type,
                        directed,
                        centered){
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # left/right objects
  Lproj <- hyp_proj$Lproj
  Rproj <- hyp_proj$Rproj
  # hypothesis set indices as a matrix
  s_ind <- matrix(FALSE,n,n)
  s_ind[hyp_indices] <- TRUE
  # projection matrix, accounting for symmetry
  if(!directed){
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind & upper.tri(s_ind,diag=TRUE)
    # linear mapping and projection basis
    Gd <- Sym_span(s_ind,s_ind_data)
    UV <- pracma::orth(crossprod(Gd,(Rproj %x% Lproj)[c(s_ind),]))
  }
  else{
    # index matrix for collecting hypothesis set from data
    s_ind_data <- s_ind
    # projection basis
    UV <- pracma::orth((Rproj %x% Lproj)[c(s_ind),])
  }
  # centering
  if(centered){
    UV <- center_orth(UV)
  }
  # final projection dimensions
  ss <- nrow(UV)
  dd <- ncol(UV)
  # test
  if(dd==1){
    X <- matrix(sapply(1:m1,function(kk){c(crossprod(UV, A[[kk]][s_ind_data]))}),ncol=1)
    Y <- matrix(sapply(1:m2,function(kk){c(crossprod(UV, B[[kk]][s_ind_data]))}),ncol=1)
  }
  else{
    X <- t(sapply(1:m1,function(kk){c(crossprod(UV, A[[kk]][s_ind_data]))}))
    Y <- t(sapply(1:m2,function(kk){c(crossprod(UV, B[[kk]][s_ind_data]))}))
  }
  # basic F test on augmented data assuming equal variance, independence
  # sse values
  sse_null <- sum((scale(rbind(X,Y),scale=F))^2)
  sse_alt <- sum((scale(X,scale=F))^2) + sum((scale(Y,scale=F))^2)
  # adjustment for tilde
  if(var_type=='quasi'){
    Xa <- (diag(ss) - tcrossprod(UV)) %*% sapply(A,function(x){x[s_ind_data]})
    Ya <- (diag(ss) - tcrossprod(UV)) %*% sapply(B,function(x){x[s_ind_data]})
    sse_alt_den <- sum(Xa^2) + sum(Ya^2)
    df_den <- (m1+m2)*(ss - dd)
  }
  else{
    sse_alt_den <- sse_alt
    df_den <- dd*(m1+m2-2)
  }
  # statistic, dfs, pvalue
  out <- list()
  out$stat <- ((sse_null - sse_alt)/dd)/(sse_alt_den/(df_den))
  out$dfs <- dd
  out$pval <- stats::pf(out$stat,df1=dd,df2=df_den,lower.tail=FALSE)
  return(out)
}


#' Mesoscale testing
#'
#' \code{Mesoscale_test} performs a hypothesis test of equality of edge
#' expectations for a prespecified collection
#' (hypothesis set) of edges between two samples of aligned networks.
#' It is compatible with either weighted (Gaussian) edge or binary edge networks
#' (specified with argument \code{edge_type}), undirected or directed networks
#' (\code{directed}), and networks with or without self loops (\code{self_loops}).
#' The mesoscale testing methodology (see MacDonald et al., 2024+) uses edge
#' variables outside the hypothesis set to learn a projection of the edges in
#' the hypothesis set, to reduce the dimension of the test and improve test power.
#'
#' @usage Mesoscale_test(A,B,sig,hyp_set,
#'                       edge_type='weighted',dimension,
#'                       directed=TRUE,self_loops=TRUE,
#'                       proj_type='impute',var_type='basic',centered,
#'                       masked_set)
#'
#' @param A a list containing networks in 1st sample; each element is an adjacency matrix.
#' @param B a list containing networks in 2nd sample, each element is an adjacency matrix, with nodes aligned with \code{A}.
#' @param sig the significance level for the hypothesis test.
#' @param hyp_set the hypothesis set of node pairs, specified as one of: (1) a 2-element list
#' containing the incident row and column indices (rectangle), (2) a list of 2-element lists each containing
#' incident row and column indices (collection of rectangles), or (3) a 2-column matrix of (row,column) index
#' pairs (unstructured hypothesis set).
#' @param edge_type a string, either \code{'weighted'} or \code{'binary'} depending on the possible
#' edge values. The test with weighted edges is performed assuming Gaussian edges. Defaults to \code{'weighted'}.
#' @param dimension a positive integer, specifies the dimension of projection in the mesoscale test.
#' @param directed a Boolean, specifies whether the network is directed. If \code{directed=FALSE}, the
#' testing methodology will be adjusted to account for symmetry. Defaults to \code{TRUE}.
#' @param self_loops a Boolean, if \code{FALSE} the test will ignore diagonal entries. Defaults to \code{TRUE}.
#' @param proj_type a string, either \code{'impute'} or \code{'one_step'} to choose the projection learning
#' approach. Recommendation is \code{'one_step'} for rectangular hypothesis sets, \code{'impute'} for unstructured
#' hypothesis sets. Defaults to \code{'impute'}.
#' @param var_type a string, either \code{'basic'} or \code{'quasi'} to choose the variance estimation approach. For
#' binary edges, \code{'quasi'} corresponds to fitting with a quasibinomial GLM. Defaults to \code{'basic'}.
#' For weighted edges, \code{'quasi'} calibrates the variance based on the variability of the edges orthogonal to the projected data, while
#' \code{'basic'} uses the variability of the projected data within the samples.
#' @param centered a Boolean, specifies whether to center the projections to ignore differences in the overall
#' edge mean on the mesoscale set. Defaults to \code{FALSE}.
#' @param masked_set a masked set of node pairs ignored in both the projection learning and testing phases
#' of the mesoscale methodology. Specified in the same format as \code{hyp_set}. Defaults to the empty set.
#'
#'
#' @return a 2-element vector consisting of the acceptance/rejection decision & p-value for the mesoscale test.
#'
#' @export
#'
#' @examples
#' A <- genSparseGraph(4,model=list(name='2SBM',n=10,p=0.5,q=0.3))
#' B <- genSparseGraph(4,model=list(name='2SBM',n=10,p=0.5,q=0.3))
#'
#' # hypothesis set specified as a rectangle
#' meso1 <- Mesoscale_test(A,B,sig=0.1,
#'                         hyp_set=list(1:4,1:4),
#'                         edge_type='binary',dimension=1,directed=FALSE,self_loops=FALSE,
#'                         proj_type='one_step',
#'                         var_type='basic')
#'
#' # same hypothesis set specified as unstructured indices
#' meso2 <- Mesoscale_test(A,B,sig=0.1,
#'                         hyp_set=cbind(c(1,1,1,2,2,3),c(2,3,4,3,4,4)),
#'                         edge_type='binary',dimension=1,directed=FALSE,self_loops=FALSE,
#'                         proj_type='one_step',
#'                         var_type='basic')
#' # check that meso1 and meso2 give the same result
#'
#' # hypothesis set and masked sets specified as rectangles
#' meso3 <- Mesoscale_test(A,B,sig=0.1,
#'                         hyp_set=list(list(1:2,1:2),list(9:10,9:10)),
#'                         edge_type='binary',dimension=1,directed=FALSE,self_loops=FALSE,
#'                         proj_type='impute',
#'                         var_type='quasi',
#'                         centered=TRUE)
#'
Mesoscale_test <- function(A,B,
                           sig,
                           hyp_set, # a rectangle (list of {row,col}), a 2-column matrix of entries, or a list if rectangles
                           edge_type='weighted',# or binary
                           # other options
                           dimension,
                           directed = TRUE,
                           self_loops = TRUE, # overrides diagonal entries in the hyp_set, masks them in the imputation routine
                           proj_type='impute',# or onestep
                           var_type='basic', # or 'quasi'
                           centered=FALSE,
                           masked_set=list(NULL,NULL) # a rectangle (list of {row,col}), a 2-column matrix of entries, or a list if rectangles
){
  # parameter checking
  # ...
  # does the test set specify a rectangle, if so expand to general indices
  if(is.matrix(hyp_set)){
    hyp_indices <- hyp_set
  }
  else{
    if(is.list(hyp_set[[1]])){
      hyp_indices <- Reduce(rbind,lapply(hyp_set,function(x){as.matrix(expand.grid(x))}))
    }
    else{
      hyp_indices <- as.matrix(expand.grid(hyp_set))
    }
  }
  # does the masked set specify a rectangle, if so expand to general indices
  if(!is.null(masked_set)){
    if(is.matrix(masked_set)){
      masked_indices <- masked_set
    }
    else{
      if(is.list(masked_set[[1]])){
        masked_indices <- Reduce(rbind,lapply(masked_set,function(x){as.matrix(expand.grid(x))}))
      }
      else{
        masked_indices <- as.matrix(expand.grid(masked_set))
        if(nrow(masked_indices)==0){
          masked_indices <- NULL
        }
      }
    }
  }
  else{
    masked_indices <- NULL
  }
  # account for self loops in masked set and hypothesis set
  if(!self_loops){
    hyp_sl <- hyp_indices[,1]==hyp_indices[,2]
    if(sum(hyp_sl) > 0){
      masked_indices <- rbind(masked_indices,hyp_indices[hyp_sl,])
      hyp_indices <- hyp_indices[!hyp_sl,]
    }
  }
  # account for symmetry in masked set and hypothesis set
  if(!directed){
    # augment hypothesis indices and masked indices to include diagonal mirrors
    hyp_indices <- unique(rbind(hyp_indices,hyp_indices[,c(2,1)]))
    if(!is.null(masked_indices)){
      masked_indices <- unique(rbind(masked_indices,masked_indices[,c(2,1)]))
    }
  }
  # dimensions
  m1 <- length(A)
  m2 <- length(B)
  n <- nrow(A[[1]])
  # group means
  Abar <- Reduce('+',A)/m1
  Bbar <- Reduce('+',B)/m2
  # transformation for binary edges
  if(edge_type=='binary'){
    delta <- 1/(3*((m1+m2)/2))
    Abar <- Logit(pmin(pmax(Abar,delta),1-delta))
    Bbar <- Logit(pmin(pmax(Bbar,delta),1-delta))
  }
  # learn left and right-hand side projections
  if(proj_type=='impute'){
    hyp_proj <- Subspace_impute(Abar,Bbar,
                                dimension,
                                hyp_indices,
                                masked_indices,
                                self_loops,
                                directed,
                                centered)
  }
  else{
    hyp_proj <- Subspace_onestep(Abar,Bbar,
                                 dimension,
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
  return(c(as.integer(test_out$pval <= sig),test_out$pval,test_out$stat))
}
