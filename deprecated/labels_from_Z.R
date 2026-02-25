# NOTE: may use this documentation later?
#' Extract Community Labels from Block Membership Matrix
#'
#' @description
#' Converts a binary block membership matrix \eqn{Z} into a vector of community
#' labels. Each row of \eqn{Z} is assumed to be a one-hot encoding of a node’s
#' community membership (i.e., exactly one entry per row is equal to 1).
#'
#' @param Z A binary membership matrix of size \eqn{n \times K}, where
#'   \eqn{n} is the number of nodes and \eqn{K} is the number of communities.
#'
#' @return An integer vector of length \eqn{n} giving the community label
#'   (from 1 to \eqn{K}) for each node.
#'
#' @examples
#' # Example: 4 nodes, 2 communities
#' Z <- matrix(c(1,0,
#'               1,0,
#'               0,1,
#'               0,1), nrow = 4, byrow = TRUE)
#' get_labels_from_Z(Z)
#' # Returns: c(1, 1, 2, 2)
#'
# get_labels_from_Z <- function(Z) {
#   apply(Z, 1, function(row) which(row == 1))
# }
