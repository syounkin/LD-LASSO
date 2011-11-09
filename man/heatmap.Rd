\name{heatmap}
\alias{heatmap}
\title{
Plots the heat map for the LD LASSO solution space
}
\description{
  This function uses the matrix beta0.mat, computed by the
  function ld_lasso_method, to plot a heat map of LD LASSO solutions.
  Increasing darkness corresponds to increasing values for the SNP LD
  LASSO estimate.
}
\usage{
heatmap( ldlasso.obj, s2.indx )
}
\arguments{
   \item{ldlasso.obj}{
    An ld lasso object
  }\item{s2.indx}{
    s2
    }
}
\details{
  See vignette.
}
\value{
  A heat map is a visual representation of the matrix of LD LASSO
  estimates.  Rows correspond to s2 values and columns correspond to SNP
  indices.  The magnitude of the SNP estimate is represented by the
  degree to which the element is shaded.  Increasing darkness
  corresponds to increasing magnitude. 
}
\author{Samuel G. Younkin}
\seealso{
  plot_ldlasso, plot_beta, image
}