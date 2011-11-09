\name{ldlasso-package}
\alias{ldlasso-package}
\alias{ldlasso}
\docType{package}
\title{LD LASSO for Case Control Genetic Association Study}
\description{
Performs a modification of the fused LASSO for SNP subset selection in
genetic association studies.  The LD LASSO requires densely spaced SNP data
and haplotype block structure.
}
\details{
The ldlasso package is a collection of functions designed to apply the LD LASSO
methodology as described in Younkin and Rao 2010.  Functions usesful for
displaying the results are also included.  
\tabular{ll}{
Package: \tab ldlasso\cr
Version: \tab 2.0\cr
Date: \tab February 1, 2011\cr
}
}
\note{
  See the vignette 'ldlasso_vignette.pdf' in ldlasso package directory
  doc/.\cr
  \cr
  For a package useful for creating haplotype block boundaries see
  MATILDE by Pattaro et al., available at http://astor.som.jhmi.edu/~gp/software/matilde/\cr
\cr
  The following function may be used to create block boudaries with MATILDE.\cr
  \cr
  find.bounds <- function( block.obj, prob.threshold = 0.95, B = 1e4, frac = 2 )\{
  Xa <- as.double.snp.data(block.obj@gtdata)\cr
  X1 <- ifelse( Xa == 0, "A", "B" )\cr
  X2 <- ifelse( Xa == 2, "B", "A" )\cr
  X <- makeGenotypes(matrix(paste(X1,X2,sep="/"), byrow = F, ncol = dim(Xa)[2] ))\cr
  LD.obj <- LD(X)\cr
  null.X <- null.density(X)\cr
  nXr <- data.frame(null.X$x, null.X$yR)\cr
  names(nXr) <- c("x", "y")\cr
  matilde.obj <- matilde(LD.obj$r^2, nXr, B=B)\cr
  n <- length(matilde.obj$L)\cr
  b.vec <- apply(matilde.obj$b[(n/frac):n,], 2, mean)\cr
  block.cood <- as.numeric((b.vec >= prob.threshold))\cr
  return( list( matilde.obj = matilde.obj, prob.threshold = prob.threshold, block.cood = block.cood ) )\cr
  \}
}
\author{
  Samuel G. Younkin\cr
  <samuel.younkin@gmail.com>
}
\references{Samuel Younkin, Joseph Nadeau, Robert Elston and J. Sunil Rao, "The Linkage Disequilibrium LASSO for SNP Selection in a Genetic Association Study of Late Onset Alzheimer Disease," Technical Report, 2010 }
\keyword{ package }
\seealso{ ld_lasso_method, ld_lasso, plot_ldlasso, plot_beta, heatmap }
\examples{
  # Load example data
  data(ldlasso_example)
  ls()

  plot_ldlasso( ldlasso.obj )

  heatmap(ldlasso.obj)

}
