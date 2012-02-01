setGeneric("findS1", function(object, ... ) standardGeneric("findS1"))

setMethod("findS1", "ldlasso", function( object, ... ) {
  if( !is.null(object@s1) )
    return(cat("s1 is not NULL\nNo changes made.\n"))
  fn.findS1(object, ... )
})

fn.findS1 <- function( ldlasso.obj, B  = 10, alpha = 0.05, tol = 5e-3, setS1 = TRUE, verbose = TRUE ){
  if(verbose) cat( paste( "Null value for s1. Finding s1 for alpha = ", alpha, "...\n", sep = "" ) )

  s1.low <- 0; s1.hi <- 10; fp.rate <- 1;

  while( abs( fp.rate - alpha ) > tol ){
    fp.tot <- 0
    s1 <- mean(c(s1.low,s1.hi))
    ldlasso.obj@s1 <- s1
    for( i in 1:B ){
      ldlasso.obj@pheno <- ldlasso.obj@pheno[sample(length(ldlasso.obj@pheno))]
      ldlasso.obj <- ldlassoSolve(ldlasso.obj)
      fp <- sum( abs( ldlasso.obj@beta ) > 1e-6 )
      fp.tot <- fp + fp.tot
    }
    fp.rate <- fp.tot/ncol(ldlasso.obj@geno)/B
    ## cat( c(fp.rate, "\n") )
    if( fp.rate < alpha ){
      ldlasso.obj@s1 <- mean(c(s1,s1.hi))
      s1.low <- s1
    }else{
      ldlasso.obj@s1 <- mean(c(s1,s1.low))
      s1.hi <- s1
    }
  }
  ldlasso.obj@s1 <- s1
  ldlasso.obj <- solve(ldlasso.obj)
  if( !setS1 ) ldlasso.obj@s1 <- NULL
  return(ldlasso.obj)
}
