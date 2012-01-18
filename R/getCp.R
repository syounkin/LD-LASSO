getCp <- function( ldlasso.obj, B = 10 ){
  logOR.obj <- logOR( geno = ldlasso.obj@geno, pheno = ldlasso.obj@pheno )
  y0 <- logOR.obj$y
  # Should add a conditional here for when beta has laready been solved
  beta0 <- solve(ldlasso.obj)
  bet0.mat <- beta0
  ystar.mat <- c()
  betastar.mat <- c()
  for( b in 1:B ){
    boot.indx <- sample(x = nrow(ldlasso.obj@geno), size = nrow(ldlasso.obj@geno), replace = TRUE)
    Xstar <- ldlasso.obj@geno[boot.indx,]
    Ystar <- ldlasso.obj@pheno[boot.indx]
    ldlasso.boot.obj <- ldlasso.obj
    ldlasso.boot.obj@geno <- Xstar
    ldlasso.boot.obj@pheno <- Ystar    
    ystar <- logOR( geno = ldlasso.boot.obj@geno, pheno = ldlasso.boot.obj@pheno )$y
    betastar <- solve( ldlasso.boot.obj )
    ystar.mat <- rbind( ystar.mat, ystar )
    betastar.mat <- rbind( betastar.mat, betastar )
  }
  a <- ( y0 - beta0 )%*%( y0 - beta0 )
  df <- 0
  for( j in 1:ncol(ystar.mat) ){
    df <- df + cov(ystar.mat[,j], betastar.mat[,j])
  }
  cp <- a - ncol(ldlasso.obj@geno) + 2*df
  return(cp)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~ Archive of old code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ( block.obj = NA, Xa = Xstar, Y = Ystar, s1 = s1, s2 = s2, r2.cut = r2.cut, form = 3, ytype = ytype, block.cood = block.cood )  
   #}else if( ytype == 2 ){
    #  Xstar = Xa
    #  Ystar = rmnorm( n = 1, mean = mean, varcov = Sigma )
    #}

    ## if( Ntot > 10 ){
    ##   if( jj%%floor(Ntot/10) == 0 ){
    ##     cat( c( 100*round(jj/Ntot,2), "% " ), sep = "" )
    ##   }
    ## }
    #if( ytype == 1 ){

  ##   Xa <- geno; Y <- pheno;
  
##   # remove SNPs in perfect LD -- to avoid noninvertible covariance matrix
##   index.mat <- which(cor(Xa)>0.99999 & !diag(ncol(Xa)), arr.ind = TRUE )
##   index.vec <- unique(index.mat[index.mat[,1]<index.mat[,2],][,2])
##   Xa <- Xa[,-index.vec]
##   y.mat <- c(); OR.mat <- c(); f0.mat <- c()
##   n.boot <- 5e3
##   for( i in 1:n.boot ){
##     boot.vec <- sample( n, size = n, replace = TRUE )
##     Xa.boot <- Xa[boot.vec,]
##     Y.boot <- Y[boot.vec]
##     logOR.obj <- logOR(Xa = Xa.boot, Y = Y.boot)
##     y.mat <- rbind( y.mat, logOR.obj$y)
##     OR.mat <- rbind( OR.mat, logOR.obj$OR )
##     f0.mat <- rbind( f0.mat, logOR.obj$f0 )
##   }
## Sigma <- var(y.mat)
## return(Sigma)

