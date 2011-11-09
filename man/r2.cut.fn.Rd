\name{r2.cut.fn}
\alias{r2.cut.fn}
\title{
Find the r-squared cutoff
}
\description{
  Returns the minimum value of r2 in r2.vec that keeps the number of constraint matrix columns less
than maxcol.
}
\usage{
r2.cut.fn(block.obj, block.cood, Xa, Y, maxcol, r2.cut.min, r2.cut.max, r2.vec.length)
}
\arguments{
  \item{block.obj}{
    An object of class gwaa.data from GenABEL.
  }
  \item{block.cood}{
  A vector of length p+1, where p is the number of SNPs.  block.cood is an indicator vector
  that indicates block boundaries at all p+1 SNP bounded intervals.  Use
  find.bounds to create this vector.
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
\item{maxcol}{
  The maximum number of columns allowed in the constraint matrix.  Default is 5,000.
}
\item{r2.cut.min}{
  The lower bound of the r2.cut interval
}
  \item{r2.cut.max}{
  The upper bound of the r2.cut interval
}
\item{r2.vec.length}{
  The number of r2.cut values uniformly spaced in r2.cut interval
}
}
\value{
  The function r2.cut.fn returns the minimum value of r2.cut that satsisfies, ncol(A) < maxcol, where A is the constraint matrix 
}
\author{
  Samuel G. Younkin
}
\seealso{
ld_lasso
}