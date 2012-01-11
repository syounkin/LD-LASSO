validldlasso = function( object ){
  if( !all( object@geno %in% 0:2 ) )
    return( "Some genotype is not equal to 0, 1 or 2" )
  if( !all( object@pheno %in% 0:1 ) )
    return( "Some phenotype is not 0 or 1" )
  if( object@s1 <= 0 )
    return( "Parameter s1 is not positive" )
  if( object@s2 <= 0 )
    return( "Parameter s2 is not positive" )
  TRUE
}

setValidity( "ldlasso", validldlasso )
