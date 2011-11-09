\name{block.cood}
\docType{data}
\alias{block.cood}
\alias{blockcood}
\title{Vector of Haplotype Block Boundaries}
\description{
  This is an indicator vector for haplotpye block boundaries.
  If there are p SNPs in the data then this vector will be of length
  p+1.  Each element in the vector indicates whether the inter-SNP
  interval is a block boundary.  The first and last element must always
  be one.
}
\usage{block.cood}
\format{A vector containing 86 entries}
\source{}
\references{Samuel Younkin, Joseph Nadeau, Robert Elston and J. Sunil
  Rao, "The Linkage Disequilibrium LASSO for SNP Selection in a Genetic
  Association Study of Late Onset Alzheimer Disease," Technical Report,
2010 }
\keyword{datasets}