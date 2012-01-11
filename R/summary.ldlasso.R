setMethod("summary", "ldlasso",
# place holder for ldlasso summary method
function( object, ... ){
  baf <- colSums(object@geno)/2/nrow(object@geno)
  sum.list <- list( baf = baf )
  cat(paste("Number of subjects: ", length(object@pheno), "\n", sep = ""))
  cat(paste("Number of SNPs: ", ncol(object@geno), "\n", sep = "" ))
  cat("\n")
  print(sum.list)
})
