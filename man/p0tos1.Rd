\name{p0tos1}
\alias{p0tos1}
\title{
Finds the LASSO parameter s1 that corresponds to desired false positive rate
}
\description{
  Estimates the LASSO parameter that corresponds to desired false
  positive rate.  More specifically, it is a bisection algorithm
  designed to find an s1 that corresponds to p0 nonzero estimates from the ld lasso with a permuted phenotype vector.
}
\usage{
p0tos1(p0, block.obj, Xa = NA, Y = NA, r2.cut = 0.01, s1high, s1low, max.iter = 100, tol = 0.1)
}
\arguments{
 \item{p0}{
    The number of nonzero estimates in permuted model.  If the expected
    false positive rate is set at 0.10, then i = 0.10*p, where p is the
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
  The value for the cutoff value of r-squared.  Can be found with r2.cut.fn. 
}
\item{s1high}{
  The initial upper limit in the bisection algorithm
}
\item{s1low}{
  The initial lower upper limit in the bisection algorithm
}
\item{max.iter}{
  The maximum number of iterations allowed in the bisection algorithm before NA is returned
}
\item{tol}{
  If p never equals i the bisection algorithim stops when |s1.old -
  s1.new| < tol.  Otherwise algortithm stops when p =i.
  }
}
\value{
  An estimate for s1
}
\author{
Samuel G. Younkin
}
\seealso{
get.s1
}