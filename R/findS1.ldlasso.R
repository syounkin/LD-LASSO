setGeneric("findS1", function(object, ... ) standardGeneric("findS1"))

setMethod("findS1", "ldlasso", function( object, ... ) {
  if( !is.null(object@s1) )
    return(cat("s1 is not NULL\nNo changes made.\n"))
  fn.findS1(object, ... )
})

fn.findS1 <- function( ldlasso.obj, loglinear = FALSE, coeff = c(-1.62, 0.75 ), iter = 10, alpha = 0.05, tol = 5e-3, setS1 = TRUE, verbose = FALSE ){
  if( !loglinear ){
    if(verbose) cat( paste( "Null value for s1. Finding s1 for alpha = ", alpha, "...\n", sep = "" ) )
  if( 2/ncol(ldlasso.obj@geno) > alpha ){
    alpha <- 1/ncol(ldlasso.obj@geno)
    cat( "decreasing alpha to 1/number of SNPs = ", alpha , "\n", sep = "" )
  }
  s1.low <- 0; s1.hi <- 5; fp.rate <- 1;
  pheno <- ldlasso.obj@pheno
  while( abs( fp.rate - alpha ) > tol ){
    fp.tot <- 0
    s1 <- mean(c(s1.low,s1.hi))
    ldlasso.obj@s1 <- s1
    for( i in 1:iter ){
      pheno.perm <- pheno[sample(length(pheno))]
      ldlasso.obj <- ldlasso(geno = ldlasso.obj@geno,
                             pheno = pheno.perm,
                             s1 = ldlasso.obj@s1,
                             s2 = ldlasso.obj@s2,
                             r2 = ldlasso.obj@r2
                             )
      ldlasso.obj <- solve(ldlasso.obj)
      fp <- sum( abs( ldlasso.obj@beta ) > 1e-6 )
      fp.tot <- fp + fp.tot
    }
    fp.rate <- fp.tot/ncol(ldlasso.obj@geno)/iter
    if(verbose) cat( c(s1, " ", fp.rate, "\n") )
    if( fp.rate < alpha ){
      ldlasso.obj@s1 <- mean(c(s1,s1.hi))
      s1.low <- s1
    }else{
      ldlasso.obj@s1 <- mean(c(s1,s1.low))
      s1.hi <- s1
    }
  }
  ldlasso.obj <- ldlasso(geno = ldlasso.obj@geno,
                             pheno = pheno,
                             s1 = s1,
                             s2 = ldlasso.obj@s2,
                             r2 = ldlasso.obj@r2
                             )
  ldlasso.obj <- solve(ldlasso.obj)
  if(!validObject(ldlasso.obj)) cat( "warning: invalid object being returned from findS1.\n")
  if( !setS1 ) ldlasso.obj@s1 <- NULL
  return(ldlasso.obj)
  }else{
    ldlasso.obj <- ldlasso(geno = ldlasso.obj@geno,
                             pheno = pheno,
                             s1 = 10^( coeff[1]+coeff[2]*log10(ldlasso.obj@s2) ),
                             s2 = ldlasso.obj@s2,
                             r2 = ldlasso.obj@r2
                             )
    ldlasso.obj <- solve(ldlasso.obj)
    if(!validObject(ldlasso.obj)) cat( "warning: invalid object being returned from findS1.\n")
    if( !setS1 ) ldlasso.obj@s1 <- NULL
    return(ldlasso.obj)
  }
}
