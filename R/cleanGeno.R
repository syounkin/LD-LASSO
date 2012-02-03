cleanGeno <- function( ldlasso.obj, ... ){
  Xa <- ldlasso.obj@geno
  mono.test <- function( geno.vec ){
    all(geno.vec==0)||all(geno.vec==2)
  }
  index.vec <- which(apply( Xa, 2, mono.test ))
#  if(length(index.vec)!=0)
#    Xa <- Xa[,-index.vec]

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
##  Y <-  ldlasso.obj@pheno
##  n <- length(Y)
##  remove SNPs in perfect LD -- to avoid noninvertible covariance matrix
## y.mat <- c(); OR.mat <- c(); f0.mat <- c()
##   n.boot <- 5e3
##   for( i in 1:n.boot ){
##     boot.vec <- sample( n, size = n, replace = TRUE )
##     Xa.boot <- Xa[boot.vec,]
##     Y.boot <- Y[boot.vec]
##     logOR.obj <- logOR(geno = Xa.boot, pheno = Y.boot)
##     y.mat <- rbind( y.mat, logOR.obj$y)
##     OR.mat <- rbind( OR.mat, logOR.obj$OR )
##     f0.mat <- rbind( f0.mat, logOR.obj$f0 )
##   }
## Sigma <- var(y.mat)
## return(Sigma)
## }
