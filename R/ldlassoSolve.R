ldlassoSolve <- function( ldlasso.obj ){

  geno <- ldlasso.obj@geno
  pheno <- ldlasso.obj@pheno
  s1 <- ldlasso.obj@s1
  s2 <- ldlasso.obj@s2
  r2 <- ldlasso.obj@r2
  delta <- 1e-6

  p <- ncol(geno)
  index.mat <- which(cor(geno)^2 > r2 & lower.tri(matrix(1, p, p)), arr.ind = TRUE )
  r <- nrow(index.mat)
  D <- matrix( 0, nrow = r, ncol = p )
  r2.vec <- c()
  for( i in 1:r ){
    D[i,index.mat[i,1]] <- -1
    D[i,index.mat[i,2]] <- 1
    r2.vec <- c( r2.vec, cor(geno[,index.mat[i,1]],geno[,index.mat[i,2]])^2 )
  }
  X2 <- geno[pheno == 1,] # 1 is case
  X1 <- geno[pheno == 0,] # 0 is control
  f2 <- colSums(X2)/2/dim(X2)[1]
  f1 <- colSums(X1)/2/dim(X1)[1]
  n0 <- sum( pheno == 0 )
  n1 <- sum( pheno == 1 )
  OR <- ifelse( ( f2 == 0 & f1 == 0) | (f2 == 1 & f1 == 1), 1, f2/(1-f2)/(f1/(1-f1)) )
  OR <- ifelse( OR == Inf, 1e6, OR )
  OR <- ifelse( OR == 0, 1/1e6, OR )
  y <- log(OR)
  var_y <- ifelse( f1 == 0 | f2 == 0, 1e6, ( n0*f1*(1-f1) + n1*f2*(1-f2) ) / ( 2*n0*n1*f1*f2*(1-f1)*(1-f2) ) )
  y <- y/sqrt(var_y)
  ldlasso.const <- -s2*log(r2.vec) + delta

  One <- rep( 1, p )
  I <- diag(1, nrow = p)
  Zero <- rep( 0, p)
  ZeroMat <- matrix( 0, nrow = p, ncol = p )
  ZeroMat_minus1 <- matrix( 0, nrow = p - 1, ncol = p )
  ZeroMat_r <- matrix( 0, nrow = r, ncol = p )
  A1 <- cbind( I, -I, I )
  A2 <- cbind( ZeroMat, I, ZeroMat )
  A3 <- cbind( ZeroMat, ZeroMat, I )
  A4 <- c( Zero, -One, -One )
  A5 <- cbind( ZeroMat_r, D, D )
  A6 <- cbind( ZeroMat_r, -D, -D )
  I_3p <- diag( 1, nrow = 3*p )
  yc = c( y, ifelse( y >= 0, y, 0 ), ifelse( y <= 0, -y, 0 ) )
  A <- t(rbind( A1, A2, A3, A4, A5, A6 ))
  b0 <- c( rep( Zero, 3 ), -s1, rep( -ldlasso.const , 2 )  )

  qp <- solve.QP(Dmat = I_3p, dvec = yc, Amat = A, bvec = b0, meq = p, factorized = FALSE)
  result <- list(  qp = qp, y = y, A = A, r2.vec = r2.vec, b0 = b0, OR = OR )

  ldlasso.obj@beta <- result$qp$solution[1:p]
  
  return(ldlasso.obj)
  
  }
