validldlasso = function( object ){
  if( !all( object@geno %in% 0:2 ) )
    return( "Some genotype is not equal to 0, 1 or 2" )
  if( !all( object@pheno %in% 0:1 ) )
    return( "Some phenotype is not 0 or 1" )
  if( !is.null(object@s1))
     if( object@s1 <= 0 )
       return( "Parameter s1 is not positive" )
  if( object@s2 <= 0 )
    return( "Parameter s2 is not positive" )
  if( length(object@beta) > 0 )
    if( (  sum(abs(object@beta)) - object@s1 ) > 1e-10 )
      return( "LASSO constraint has been violated." )
  ## if( length(object@beta) > 0 ){
  ##     p <- ncol(object@geno)
  ##     index.mat <- which(cor(object@geno)^2 > object@r2 & lower.tri(matrix(1, p, p)), arr.ind = TRUE )
  ##     if( any( abs(abs(object@beta[index.mat[,1]]) - abs(object@beta[index.mat[,2]]) ) - ( -object@s2*log(cor(object@geno)[index.mat]^2) + 1e-6 ) > 1e-10 ) )
  ##     return ("LD LASSO constraint violated")
  ##   }
  TRUE
}

setValidity( "ldlasso", validldlasso )
