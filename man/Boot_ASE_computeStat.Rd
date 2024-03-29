% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Boot_ASE_EPA_test.R
\name{Boot_ASE_computeStat}
\alias{Boot_ASE_computeStat}
\title{Boot-ASE/EPA Compute Statistics}
\usage{
Boot_ASE_computeStat(A_set, B_set, r)
}
\arguments{
\item{A_set}{a sparse EPAacency matrix}

\item{B_set}{a sparse EPAacency matrix}

\item{r}{scalar specifying rank of population EPAacency}
}
\value{
ASE statistic and EPA statistic
}
\description{
This function computes the test statistics for Boot-ASE/EPA test, serving as an auxiliary functions for the function: LowRankTests
}
\examples{
A <- genSparseGraph(4,model=list(name='2SBM',n=4,p=0.3,q=0.3))
B <- genSparseGraph(4,model=list(name='2SBM',n=4,p=1,q=0.3))
Boot_ASE_computeStat(A, B, 2)

}
