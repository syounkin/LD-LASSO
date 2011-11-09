\name{plot_beta}
\alias{plot_beta}
\title{
Plot the LD LASSO solution by physical position.
}
\description{
  Creates a plot of LD LASSO estimates using the physical map, bpmap
}
\usage{
plot_beta( beta, bpmap=NA, block.bounds.vec=NA )
}
\arguments{
   \item{beta}{
    A vector of ldlasso estimates.
  }
  \item{bpmap}{
    A vector of map positions in terms of kb from the left boundary
  }
  \item{block.bounds.vec}{
    A vector of block boundaries positions in kilobases from the left boundary.
    Boundaries are plotted at the midpoint of SNP intervals.
  }
  }
\details{
 Includes haplotpye block boundaries as vertical broken lines.
}
\value{
  A plot object.
}
\author{Samuel G. Younkin}
\seealso{
  plot_ldlasso
}