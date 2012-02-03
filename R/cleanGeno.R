cleanGeno <- function( ldlasso.obj, ... ){
  Xa <- ldlasso.obj@geno
  mono.test <- function( geno.vec ){
    all(geno.vec==0)||all(geno.vec==2)
  }
  index.vec <- which(apply( Xa, 2, mono.test ))
  index.mat <- which(cor(Xa)>0.99999 & !diag(ncol(Xa)), arr.ind = TRUE )
  if(!(nrow(index.mat)==0)){
    index.mat <- index.mat[diff(index.mat)>0,]
    if(is.null(dim(index.mat))){
      index.vec <- c( index.vec, index.mat[2] )
    }else{
      index.vec <- c( index.vec, index.mat[,2] )
    }
    Xa <- Xa[,-index.vec]
  }
  if( length(index.vec)!=0 ){
    cat( "Removed SNPs with indices", index.vec, "\n", sep = "" )
    ldlasso.obj <- ldlasso(geno = Xa,
                         pheno = ldlasso.obj@pheno,
                         s1 = ldlasso.obj@s1,
                         s2 = ldlasso.obj@s2,
                         r2 = ldlasso.obj@r2
                         )
    return(ldlasso.obj)
  }else{
    cat("Geno is already clean.\n")
    return(ldlasso.obj)
  }
}
