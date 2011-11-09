\name{block.map.matrix}
\alias{block.map.matrix}
\title{ Creates an indicator matrix for haplotype block boundaries, for
  use in ld_lasso.
}
\description{
Simple function that maps the block boundary vector to an indicator
matrix for use in the definition of constraint matrix.  This matrix
ensures that only within block SNP pairs are considered. 
}
\usage{
block.map.matrix(block.cood)
}
\arguments{
  \item{block.cood}{
  A vector of length p+1, where p is the number of SNPs.  block.cood is an indicator vector
  that indicates block boundaries at all p+1 SNP bounded intervals.  Use
  find.bounds to create this vector.
  }
}
\value{
  A matrix of logical variables.  If the (i,j) entry is TRUE than SNP i
  and SNP j are in the same haplotype block.
}

\author{
  Samuel G. Younkin
}
\seealso{
ld_lasso
}
