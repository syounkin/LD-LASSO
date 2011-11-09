\name{ldlasso.obj}
\docType{data}
\alias{ldlasso.obj}
\alias{ldlassoobj}
\title{LD LASSO Object created by the function ld_lasso_method}
\description{ This is a data frame created by the function
  ld_lasso_method.  It contains information about the parameters used to
  generate the three LD LASSO solutions, as well as the solutions
  themselves.  It also contains the results from the grid search
  performed to estimates the cp-optimal solution. (see cp.obj)  Some
  basic map information is also stored here for use when plotting.
}
\usage{ldlasso.obj}
\format{data.frame\cr
  \tabular{lll}{
    beta1 \tab numeric \tab cp-optimal solution \cr
    beta2 \tab numeric \tab fused solution \cr
    beta3 \tab numeric \tab un-fused solution \cr
    s2star \tab numeric \tab s2 parameter used for cp-optimal solution \cr
    cp.obj \tab data.frame \tab s2 and cp vectors along with solution
    matrix used to find s2star \cr
    log10p \tab numeric \tab log base 10 of p-values from
    allelic test of association for each SNP \cr
    bpmap \tab numeric \tab vector of physical SNP positions in kb from
    left boundary \cr
    block.bounds.vec \tab numeric \tab position of block boundaries in
    kb from left boundary (does not inlcude endpoints) \cr
    s1 \tab numeric \tab LASSO parameter \cr
    B \tab numeric \tab number of bootstrap samples used by
    ld_lasso_method to find s2star \cr
    s2.vec.length \tab numeric \tab length of vector of s2 values used
    in ld_lasso_method \cr
    }
  }
\source{}
\references{Samuel Younkin, Joseph Nadeau, Robert Elston and J. Sunil
  Rao, "The Linkage Disequilibrium LASSO for SNP Selection in a Genetic
  Association Study of Late Onset Alzheimer Disease," Technical Report,
2010 }
\keyword{datasets}