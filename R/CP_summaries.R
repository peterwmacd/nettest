#### helper functions ####

# compute scan statistics from k-hop neighbourhoods
# G is an igraph object
scan_stat <- function(G, k) {
  max(sapply(igraph::V(G), function(v) {
    nbhd <- igraph::neighborhood(G, order = k, nodes = v)[[1]]
    subG <- igraph::induced_subgraph(G, nbhd)
    igraph::ecount(subG)
  }))
}

# compute 9 summary features based a adjacency matrix A, binary
# undirected with no self-loops (diagonals all 0)
summary_features <- function(A){
  # dimension
  n <- nrow(A)
  # data cleaning, G is an igraph object
  G <- igraph::graph_from_adjacency_matrix(A,mode='undirected')
  # allocate space
  fvec <- numeric(9)

  # f1: total number of edges 'size'
  # normalized by number of potential edges
  fvec[1] <- sum(pracma::triu(A,1)) / choose(n,2)

  # f2: maximum node degree 'maxdeg'
  fvec[2] <- max(rowSums(A))

  # f3: largest eigenvalue 'maxeig'
  fvec[3] <- irlba::partial_eigen(A,1)$values

  # f4-f6: 1,2,3 hop scan statistics 'scan-1','scan-2','scan-3'
  # normalized by number of potential edges
  fvec[4:6] <- sapply(1:3,function(kk){scan_stat(G,kk)}) / choose(n,2)

  # f7: triangle count 'triangle'
  # normalized by number of potential triangles
  fvec[7] <- sum(igraph::count_triangles(G)) / choose(n,3)

  # f8: global clustering coef 'gclust'
  fvec[8] <- igraph::transitivity(G,type='global')

  # f9: negated average path length 'apath'
  if (igraph::is_connected(G)) {
    fvec[9] <- -igraph::mean_distance(G)
  } else {
    dists <- igraph::distances(G)
    n <- igraph::vcount(G)
    max_dist <- max(dists[is.finite(dists)])
    dists[!is.finite(dists)] <- 2 * max_dist
    fvec[9] <- -mean(dists[lower.tri(dists)])
  }

  # return result
  return(fvec)
}

#### main function ####

#' Graph Summaries (Invariants) for Changepoint Detection
#'
#' Computes 9 graph summary statistics (or graph invariants), described in \href{https://ieeexplore.ieee.org/document/6380528}{Park et al., (2012)}.
#' for an ordered sequence of undirected networks with no self-loops.
#' 1. \code{density}: Edge density, number of edges divided by number of node pairs.
#' 2. \code{maxdeg}: Maximum node degree.
#' 3. \code{maxeig}: Maximum eigenvalue of adjacency matrix.
#' 4. \code{scan-1}: Maximum number of edges over 1-hop neighbourhoods.
#' 5. \code{scan-2}: Maximum number of edges over 2-hop neighbourhoods.
#' 6. \code{scan-3}: Maximum number of edges over 3-hop neighbourhoods.
#' 7. \code{triangle}: Triangle density, number of triangles divided by number of node triples.
#' 8. \code{gclust}: Global clustering coefficient, a measurement of triangle closure.
#' 9. \code{apath}: Negative of average path length over all node pairs. Disconnected node pairs are assigned a path length equal to twice the maximum observed (finite) path lengths.
#'
#' @param A A list of matrices, ordered sequence of adjacency matrices.
#'
#' @return A \code{length(A)}\eqn{\times}\code{9} matrix with named columns, with the 9 summary statistics computed for each adjacency matrix.
#'
#' @export
#'
#' @examples
#' CP <- list(type='kidneyegg',t_start=5,n_egg=20,p=0.5,q=0.8)
#' A <- Simulate_CP(100,10,CP=CP)$A
#' test <- CP_summaries(A)
CP_summaries <- function(A){
  # initial cleaning and dimensions, make A,B one element lists if they are not
  A <- checklist(A)
  m <- length(A)
  # compute summaries
  summaries <- t(apply(list_to_array(A),3,summary_features))
  colnames(summaries) <- c('density','maxdeg','maxeig',
                           'scan-1','scan-2','scan-3',
                           'triangle','gclust','apath')
  rownames(summaries) <- 1:m
  return(summaries)
}
