\name{ld_lasso}
\alias{ld_lasso}
\title{The main function of the LD LASSO method}
\description{
  The LD LASSO uses the correlation of SNP genotypes in a penalized least squares regression framework.  The estimator
  is the solution to a convex optimization problem, and here we use the
  solution from the package quadprog.
}
\usage{
ld_lasso(block.obj, block.cood = NA, Xa = NA, Y = NA, s1, s2, r2.cut = 0.5, delta =
1e-10, form = 3, ytype = 1, solve = TRUE )
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
\item{s1}{
  The LASSO constraint parameter -- the sum of the magnitude of the
  estimates is bounded by s1.
}
\item{s2}{
  The LD LASSO constraint parameter -- the absolute difference of SNP pair
  estimates is bounded by s2 times the log of r-squared
}
\item{r2.cut}{
  Only SNP pairs with correlation greater than r2.cut are bounded by the
  LD LASSO constraint.

}
\item{delta}{
  Included so that optimization is numerically feasible in cases when
  r-squared = 1
}
\item{form}{
  Form of constraint matrix. form is either 1, 2 or 3:\cr
  1 for cpcc.vec <- 1e6*rep(1,length(r2)) -- LASSO solution\cr
  2 for cpcc.vec <- -s2*log(r2) + delta, s1 <- 1e6 -- LD fused solution\cr
  3 for cpcc.vec <- -s2*log(r2) + delta -- LD LASSO
}
\item{ytype}{
  If ytype is 1 then Y is a vector of binary disease phenotypes, 0 for
  non-disease, 1 for diseased.  If ytype is 2 then Y is the normalized log OR.
}
\item{solve}{
  logical variable indicating whether or not to solve regression
  problem.  Useful when ld_lasso is used to construct constraint matrix,
  and the solution is not necessary, as in the selction of the r2 cutoff. 
}
}
\details{
  This function performs the ld lasso regression with parameters s1, s2
  and r2.cut on block.obj with haplotype block boundaries defined
  by block.cood.
}
\value{
  \item{qp}{
    List from the function solve.QP in the package quadprog.  This object
    contains the solutions for c( beta, beta+, beta- ) and so
    the LD LASSO estimates are the first p elements, or qp$solution[1:p]
  }
  \item{y}{
    A vector of normalized log OR
  }
  \item{A}{
    The constraint matrix
  }
  \item{r2}{
    The vector of r-squared values used to define constraint matrix.
    Elements in this vector are the correlation estimates for all
    inter-block SNP pairs.
  }
  \item{b0}{
    The vector of ld lasso constraints with length 3p
  }
  \item{OR}{
    A vector of odds ratios
  }  
}
\references{
D. Goldfarb and A. Idnani, "A numerically stable dual method for solving strictly
convex quadratic programs," Mathematical Programming, vol. 27, pp. 1-33, 1983.
}
\author{
Samuel G. Younkin
}
\seealso{
ld_lasso_method, quadprog, solve.QP
}
\examples{
  data("ldlasso_example")
  ldlasso.test <- ld_lasso( block.obj, block.cood, s1 = 1, s2 = 0.5 )

}

