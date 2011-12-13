setGeneric("geno", function(object) standardGeneric("geno"))
setMethod("geno", signature("ldlasso"), function(object) object@geno )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setGeneric("pheno", function(object) standardGeneric("pheno"))
setMethod("pheno", signature("ldlasso"), function(object) object@pheno )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setGeneric("parameters", function(object) standardGeneric("parameters"))
setMethod("parameters", signature("ldlasso"), function(object){
  return(list( s1 = object@s1, s2 = object@s2, r2 = object@r2, delta = object@delta ))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setGeneric("indexmat", function(object) standardGeneric("indexmat"))
setMethod("indexmat", signature("ldlasso"), function(object){
    return(which(cor(geno(object))^2 > parameters(object)$r2 & lower.tri(matrix(1, ncol(geno(object)), ncol(geno(object)))), arr.ind = TRUE ))
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setGeneric("beta", function(object) standardGeneric("beta"))
setMethod("beta", signature("ldlassoSolution"), function(object) object@beta )
