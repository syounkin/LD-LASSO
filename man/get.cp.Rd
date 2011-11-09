\name{get.cp}
\alias{get.cp}
\title{
Uses bootstrap sampling over a vector of LD LASSO constraint parameters,
s2, to compute a vector of cp estimates.
}
\description{
The vector of cp estimates is used to identify the cp-optimal solution.
}
\usage{
get.cp(s2low, s2high, s2.vec.length, block.obj, Xa = NA, Y = NA, s1, r2.cut, block.cood, B = 20)
}
\arguments{
\item{s2low}{
  The lower limit for the s2 vector.
}
\item{s2high}{
  The upper limit for the s2 vector.
}
\item{s2.vec.length}{
  The number of exponentially spaced values in the s2 vector.
}
\item{block.obj}{
    An object of class gwaa.data from GenABEL.
}
\item{Xa}{
    If block.obj is NA then a genotype matrix must be provided.  Xa is a
    matrix of genotype values codes as 0, 1 or 2 for homozygous major,
    heterozygous, or homozygous minor, respectively.
}
\item{Y}{
    If block.obj is NA then a phenotype vector Y must be provided.  Y is a
    vector of diagnoses, where 0 is non-diseased and 1 is diseased.
}
\item{s1}{
  The LASSO parameter
}
\item{r2.cut}{
  Only SNP pairs with correlation greater than r2.cut are bounded by the
  LD LASSO constraint.
}
\item{block.cood}{
  A vector of length p+1, where p is the number of SNPs.  block.cood is an indicator vector
  that indicates block boundaries at all p+1 SNP bounded intervals.  Use
  find.bounds to create this vector.
}
\item{B}{
    Number of bootstrap samples
}
}
\value{
  \item{s2.vec}{
    A vector of s2 values
  }
  \item{cp.vec}{
    A vector of cp estimates
  }
  \item{beta0.mat}{
    A matrix of LD LASSO estimates
  }
  \item{s1}{
    The LASSO parameter
  }
}
\author{
  Samuel G. Younkin
}

\seealso{
ld_lasso_method
}