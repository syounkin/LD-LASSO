getSigma <- function( ldlasso.obj ){

  Xa <- ldlasso.obj@geno
  Y <-  ldlasso.obj@pheno
  n <- length(Y)
  # remove SNPs in perfect LD -- to avoid noninvertible covariance matrix
  index.mat <- which(cor(Xa)>0.99999 & !diag(ncol(Xa)), arr.ind = TRUE )
  index.vec <- unique(index.mat[index.mat[,1]<index.mat[,2],][,2])
  if( length(index.vec) > 0 )
    Xa <- Xa[,-index.vec]
  y.mat <- c(); OR.mat <- c(); f0.mat <- c()
  n.boot <- 5e3
  for( i in 1:n.boot ){
    boot.vec <- sample( n, size = n, replace = TRUE )
    Xa.boot <- Xa[boot.vec,]
    Y.boot <- Y[boot.vec]
    logOR.obj <- logOR(geno = Xa.boot, pheno = Y.boot)
    y.mat <- rbind( y.mat, logOR.obj$y)
    OR.mat <- rbind( OR.mat, logOR.obj$OR )
    f0.mat <- rbind( f0.mat, logOR.obj$f0 )
  }
Sigma <- var(y.mat)
return(Sigma)
}
