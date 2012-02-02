logOR <- function( geno, pheno){
  Xa <- geno; Y <- pheno;
  X2 <- Xa[Y == 1,] # 1 is case
  X1 <- Xa[Y == 0,] # 0 is control
  f2 <- colSums(X2)/2/dim(X2)[1]
  f1 <- colSums(X1)/2/dim(X1)[1]
  f0 <- colSums(Xa)/2/dim(Xa)[1]
  n0 <- sum( Y == 0 )
  n1 <- sum( Y == 1 )
  OR <- ifelse( ( f2 == 0 & f1 == 0) | (f2 == 1 & f1 == 1), 1, f2/(1-f2)/(f1/(1-f1)) )
  OR <- ifelse( OR == Inf, 1e6, OR )
  OR <- ifelse( OR == 0, 1/1e6, OR )
  y <- log(OR)
  var_y <- ifelse( f1 == 0 | f2 == 0, 1e6, ( n0*f1*(1-f1) + n1*f2*(1-f2) ) / ( 2*n0*n1*f1*f2*(1-f1)*(1-f2) ) )
  y <- y/sqrt(var_y)
  return(list( OR = OR,
              y = y,
              var_y = var_y,
              f1 = f1,
              f2 = f2,
              f0 = f0,
              X2 = X2,
              X1 = X1,
              n0 = n0,
              n1 = n1 )
         )
}
