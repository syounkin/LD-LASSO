\name{trace.plot}
\alias{trace.plot}
\title{
Creates a trace plot for the LD LASSO
}
\description{This function creates a trace plot as a function of s2, the
  LD LASSO constraint
}
\usage{
trace.plot(beta0.mat, s2.vec, type = "l", indx = NA, s2star = NA, abs = TRUE)
}
\arguments{
  \item{beta0.mat}{
    A matrix where rows are ld lasso estimates.  This matrix has number of rows equal
    to length(s2.vec) and number of columns equal to the number of SNPs
}
\item{s2.vec}{
  The vector of s2 values used to generate the rows of beta0.mat
}
\item{type}{
  Type of plot, "l" for lines, "p" for points
}
\item{indx}{
  If index is NA then indx, is a vector of SNP indices to be
  plotted in red.
}
\item{s2star}{
  Value of s2star from cp object.  If s2star is not NA, a vertical line is plotted at this value.
}
\item{abs}{
  A logical variable indicating whether to plot absolute value of beta instead of beta
}
}
\details{
  In this trace plot each line represents a SNP estimate and we see the
  behavior of the solution as the LD LASSO constraint is varied. 
}
\value{
  Creates a trace plot.
}
\author{
  Samuel G. Younkin
}