# code to generate from different (independent edge) graph models

library(igraph)

#' @export
genSparseGraph <- function(m, model) {
  # m: the number of graphs
  A <- list()
  for (i in 1:m) {
    switch(model$name,
           ER = {
             # A1 <- triu(ceiling(sprand(model$n, model$n, model$p)), 1)
             A1 <- pracma::triu(matrix(stats::rbinom(model$n^2, 1, model$p), model$n, model$n),k=1)
             A[[i]] <- A1 + t(A1)
           },
           `2SBM` = {
             n1 <- floor(model$n/2)
             n2 <- model$n - n1
             A11 <- pracma::triu(matrix(stats::rbinom(n1^2, 1, model$p), n1, n1),k=1)
             A22 <- pracma::triu(matrix(stats::rbinom(n2^2, 1, model$p), n2, n2),k=1)
             A12 <- matrix(stats::rbinom(n1*n2, 1, model$q), n1, n2)
             A1 <- rbind(cbind(A11, A12), cbind(matrix(0, n2, n1), A22))
             A[[i]] <- A1 + t(A1)
           },
           SBM = {
             inc <- floor(model$n / model$k)
             A1 <- matrix(0, nrow = model$n, ncol = model$n)
             ni <- 0
             for (ci in 1:model$k) {
               nj <- ni
               A1[(ni + 1):(ni + inc), (nj + 1):(nj + inc)] <- pracma::triu(matrix(stats::rbinom(inc^2, 1, model$p), inc, inc),k=1)
               for (cj in (ci + 1):model$k) {
                 nj <- nj + inc
                 A1[(ni + 1):(ni + inc), (
                   nj + 1):(nj + inc)] <- pracma::triu(matrix(stats::rbinom(inc^2, 1, model$q), inc, inc),k=1)
               }
               ni <- ni + inc
             }
             A[[i]] <- A1 + t(A1)
           },
           IER = {
             EA <- pracma::triu(model$P,k=1)
             n = model$n
             A1 <- matrix(stats::runif(n*n) < EA, n, n)*1
             A[[i]] <- A1 + t(A1)
           },
           {
             stop("Model name unavailable.")
           })
  }
  return(A)
}
