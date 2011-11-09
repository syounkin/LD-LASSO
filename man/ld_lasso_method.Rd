\name{ld_lasso_method}
\alias{ld_lasso_method}
\title{
The LD LASSO function
}
\description{
  This function implements a method for the automatic selection of
  parameters for the LD LASSO.  It returns three
  solutions, the fused, cp-optimal, and unfused solutions.  It also
  creates a matrix of solutions needed for
  creating the trace plot. 
}
\usage{
ld_lasso_method(block.obj, block.cood = NA, Xa = NA, Y = NA, bpmap = NA,
maxcol = 5e3, p.frac = 0.10, B = 5, s2low = 5e-3, s2high = 5e1, s2.vec.length = 4, null = FALSE)
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
  \item{bpmap}{
    A vector of map positions in terms of kB from the left boundary
    }
\item{maxcol}{
  The upper limit on the number of
  columns in the constraint matrix.  This is to prevent computational
  overload.  Increasing maxcol may increase computation time and
  memory needed. 
}
\item{p.frac}{
  The fraction of SNPs allowed in LASSO model under null hypothesis.  This
  parameter is used in the function get.s1.    
}
\item{B}{
  The number of bootstrap iterations for cp estimate.
}
\item{s2low}{
  The lower limit of the s2 vector
}
\item{s2high}{
  The upper limit of the s2 vector
}
\item{s2.vec.length}{
  The number of exponentially spaced values of s2.
}
\item{null}{
  A logical variable that indicates if analysis should be performed on
  permuted phenotype vector.
}
}
\details{
  Use function find.bounds to create block.cood vector with MATILDE MCMC
  methods. See ldlasso help page for details on the package MATILDE.
}
\value{
    \item{beta1}{
    The ld lasso solution with s2 that minimizes cp (cp-optimal solution)
  }
  \item{beta2}{
    The ld lasso solution for lower limit of s2 interval (fused
    solution)
  }\item{beta3}{
    The ld lasso solution for upper limit of s2 interval (unfused solution)
  }
  \item{s2star}{
    The s2 value that minimizes cp
  }
  \item{cp.obj}{
    A list that contains information used for cp estimation
  }
  \item{log10p}{
    A vector of log10 p values for test of allelic association
  }
  \item{bpmap}{
    A vector of map positions in base pairs
  }
  \item{block.bounds.vec}{
    A vector of block boundaries in kB from left boundary
  }
  \item{s1}{
    the LASSO constraint
  }
  \item{B}{
    Number of bootstrap samples
  }
  \item{s2.vec.length}{
    The length of the s2 vector
    }
}
\references{Samuel Younkin, Joseph Nadeau, Robert Elston and J. Sunil Rao, "The Linkage Disequilibrium LASSO for SNP Selection in a Genetic Association Study of Late Onset Alzheimer Disease," Technical Report, 2010 }
\author{Samuel G. Younkin}
\seealso{
  ld_lasso
}
\examples{
  # Load example data
  # data(ldlasso_example)
  
  # Run the method with low B and s2.vec.length first to test.
  # ldlasso.obj <- ld_lasso_method( block.obj, block.cood, B = 3, s2.vec.length = 2 )

}