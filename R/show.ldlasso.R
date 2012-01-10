setMethod("show", "ldlasso",
function(object){
  # a place holder for the show method of ldlasso
  cat( "s1 =", object@s1, "\n" )
  cat( "s2 =", object@s2, "\n" )
  cat( "r2 =", object@r2, "\n" )
  cat( "\n" )
  cat("genotypes:\n")
  print(head(object@geno))
  if(nrow(object@geno) > 6)
    cat("...\n\n")
  cat("phenotypes:\n")
  print(object@pheno)
  cat("\n")
  cat("beta:\n")
  print(object@beta)
})
