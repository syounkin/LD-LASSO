\name{chi2}
\alias{chi2}
\title{
Creates a vector of chi-squared values from the allelic test of association
}
\description{
One degree of freedom Chi-squared test of allelic counts
}
\usage{
chi2(Xa, Y)
}
\arguments{
  \item{Xa}{
    A matrix of genotype values  coded as 0, 1 or 2 for homozygous major,
    heterozygous, or homozygous minor, respectively.  Rows are subjects and
    columns are SNPs.  May be created with method as.double.snp.data from the GenABEL
    package.
 }
\item{Y}{
  A vector of disease phenotypes, diseased = 1, non-diseased = 0
}
}
\value{
  This function returns a vector of chi-squared, one degree of freedom,
  statistics for the allelic test of association.  The length of this
  vector is the number of SNPs.
}
\author{
  Samuel G. Younkin
}
