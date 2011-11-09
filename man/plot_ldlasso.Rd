\name{plot_ldlasso}
\alias{plot_ldlasso}
\alias{plot_ld_lasso}
\alias{plot_ld_lasso_method}
\title{
Diagnostic plot for an ldlasso object
}
\description{
Creates four plots; Top-left: Allelic test of association; Top-
right: cp-optimal, fused, and unfused solutions (black, red, green); Bottom-left:
Trace plot with fixed s1 ; Bottom-right: cp estimates vs. s2
}
\usage{
plot_ldlasso( ldlasso.obj )
}
\arguments{
   \item{ldlasso.obj}{
    An ld lasso object
}
}
\details{
  Note that often two solutions overlap, and one is therefore concealed.  The
  vertical line in the cp plot is drawn at the value of s2 used in the
  cp-optimal solution.
}
\value{
  A plot object
}
\author{Samuel G. Younkin}
\seealso{
  ld_lasso_method, ld_lasso
}