% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Boot_spectral_test.R
\name{ShufflingTests}
\alias{ShufflingTests}
\title{Shuffling Tests}
\usage{
ShufflingTests(A, B, sig, bs)
}
\arguments{
\item{A}{cell array containing network in 1st population; each cell is a sparse adjacency matrix}

\item{B}{cell array containing network in 2nd population, each cell being a sparse adjacency matrix}

\item{sig}{significance level for acceptance of null hypothesis}

\item{bs}{bootstrap times}
}
\value{
tFro: output of Boot-Frobenius test (1 if null is rejected, 0 otherwise)\cr
tOp: output of Boot-Spectral test (1 if null is rejected, 0 otherwise)\cr
pvalFro: p-value for Boot-Frobenius test\cr
pvaOP: p-value for Boot-Spectral test
}
\description{
This function performs the Boot-Frobenius and Boot-Spectral tests, which are two-sample testing procedures for large graphs.
The Boot-Frobenius test is a two-sample testing procedure for large graphs. The test is based on the Frobenius norm and uses bootstrapping to approximate the null distribution of a test statistic.
The Boot-Spectra is based on adjacency spectral embedding (ASE) and uses bootstrapping to approximate the null distribution of a test statistic. \cr\cr
The test involves computing the rejection decision and the p-value, which is computed to ensure that the null hypothesis (i.e., both populations have the same population adjacency) is rejected only when the test statistic is in the upper α-quantile for both approximate distributions.
}
\examples{
A <- genSparseGraph(4, model=list(name='ER',n=10,p=.1))
B <- genSparseGraph(4, model=list(name='ER',n=10,p=.5))
test_result = ShufflingTests(A, B, 0.05, 200)

}
