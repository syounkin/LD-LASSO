setMethod("show", "ldlasso",
function(object){
  # a place holder for the show method of ldlasso
  cat( "s1 = ", object@s1, ", ", sep = "" )
  cat( "s2 = ", object@s2, ", ", sep = "" )
  cat( "r2 = ", object@r2, "\n", sep = "" )
  cat( "\n" )
  cat( nrow(object@geno), "subjects and", ncol(object@geno), "SNPs\n" )
  if( length(object@beta) == 0 ){
    cat( "\nNo solution.\n")
  }else{
    cat("\nSummary of solution vector (beta):\n")
    print(summary(object@beta))
  }
  ## cat("genotypes:\n")
  ## print(head(object@geno))
  ## if(nrow(object@geno) > 6)
  ##   cat("...\n\n")
  ## cat("phenotypes:\n")
  ## print(object@pheno)
  ## cat("\n")
  ## cat("beta:\n")
  ## print(object@beta)
})
