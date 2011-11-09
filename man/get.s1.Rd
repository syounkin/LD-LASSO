\name{get.s1}
\alias{get.s1}
\title{
Finds an estimate for s1 based on an estimate of the expected false positive
rate.  This function is a wrapper for p0tos1.
}
\description{
  This function calls the function p0tos1, which finds an estimate for
  s1, given p0, where p0 is the number of nonzero estimates found with a permuted
  phenotype vector.
}
\usage{
get.s1(i, block.obj, Xa = NA, Y = NA, r2.cut, s1high, s1low, max.iter, tol)
}
\arguments{
  \item{i}{
    Number of nonzero estimates in permuted model.  If the expected
    false positive rate is desired to be 0.10, then set i = 0.10*p, where p is the
    number of SNPs.
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
  \item{r2.cut}{
  Value for the cutoff value of r-squared.  Can be found with r2.cut.fn. 
}
\item{s1high}{
  The initial upper limit in the bisection algorithm
}
\item{s1low}{
  The lower upper limit in the bisection algorithm
}
\item{max.iter}{
  The maximum number of iterations allowed in the bisection algorithm before NA is returned
}
\item{tol}{
  If p never equals i the bisection algorithim stops when |s1.old - s1.new| < tol.  Otherwise algortithm stops when p = i.
}
}
\details{This function increases the value for p0 if NA is returned from
  p0tos1.  The bisection algorthim is contained in the function p0tos1
}
\value{
  Given the number of nonzero SNPs allowed with permuted phenotype
  vector, p0, this function returns an estimate for s1.
}
\author{
Samuel G. Younkin
}
\seealso{
  p0tos1
}