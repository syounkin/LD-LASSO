setClass("ldlasso", representation = representation( geno = "matrix", pheno = "numeric", s1 = "numeric", s2 = "numeric", r2 = "numeric", delta = "numeric" ), validity = function( object ) {
  if( object@s1 < 0 )
    return("s1 is less than zero")
  if( object@s2 < 0 )
    return("s2 is less than zero")
  if( object@r2 < 0 )
    return("r2 is less than zero")
  if( object@delta < 0 | object@delta > 1e-6)
    return("delta is either negative or too large (>1e-6)")
  if( !all(geno(object) %in% 0:2 ) )
    return("Some genotype is not in (0, 1, 2)")
  if( !all(pheno(object) %in% 0:1 ) )
    return("Some phenotype is not in (0, 1)")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setClass("ldlassoSolution", contains = "ldlasso", representation = representation( beta = "numeric" ) )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if( any( abs(abs(beta(object)[indexmat(object)[,1]]) - abs(beta(object)[indexmat(object)[,2]]) ) > ( -parameters(object)$s2*log(cor(geno(object))^2)[indexmat(object)] + parameters(object)$delta ))) # DON'T FORGET ABOUT DELTA
#	if( sum(abs(object@beta)) > (object@s1+1e-8) )
#          return("LASSO constraint violated")
## setValidity("ldlassoSolution", function(object){
##   
##           return("LD LASSO constraint violated")
## })
# LD lASSO constraint
#if( any( abs(abs(object@beta[object@indexmat[,1]]) - abs(object@beta[object@indexmat[,2]]) ) > ( -object@s2*log(cor(object@geno)[object@indexmat]^2) + object@delta ) ) )

