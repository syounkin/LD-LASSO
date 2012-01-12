ldlassoModel <- function( ldlasso.obj ){
  if( is.null(ldlasso.obj@s1) ) {
    ldlasso.obj@s1 <- findS1( ldlasso.obj )
  }
  ldlasso.obj@beta <- solve(ldlasso.obj) 
  return( ldlasso.obj)
}

findS1 <- function( ldlasso.obj, B = 1e1, alpha = 0.05, tol = 1e-3 ){
  cat( paste( "Null value for s1. Finding s1 for alpha = ", alpha, "...\n\n", sep = "" ) )
  s1.low <- 0; s1.hi <- 10; fp.rate <- 1;  
  while( abs( fp.rate - alpha ) > tol ){
    fp.tot <- 0
    s1 <- mean(c(s1.low,s1.hi))
    ldlasso.obj@s1 <- s1
    for( i in 1:B ){
      ldlasso.obj@pheno <- ldlasso.obj@pheno[sample(length(ldlasso.obj@pheno))]
      ldlasso.obj@beta <- solve(ldlasso.obj)
      fp <- sum( abs( ldlasso.obj@beta ) > 1e-6 )
      fp.tot <- fp + fp.tot
    }
    fp.rate <- fp.tot/ncol(ldlasso.obj@geno)/B
    cat( c(fp.rate, "\n") )
    if( fp.rate < alpha ){
      ldlasso.obj@s1 <- mean(c(s1,s1.hi))
      s1.low <- s1
    }else{
      ldlasso.obj@s1 <- mean(c(s1,s1.low))
      s1.hi <- s1
    }
  }
  return(ldlasso.obj@s1)
}
