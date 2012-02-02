ldlassoSolve <- function( ldlasso.obj ){

  if( is.null(ldlasso.obj@s1) )
    return( "Cannot solve!  The parameter s1 is NULL" )

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

  X2 <- ldlasso.obj@X.case
  X1 <- ldlasso.obj@X.con
  f2 <- ldlasso.obj@maf.case
  f1 <- ldlasso.obj@maf.con
  n0 <- ldlasso.obj@n.con
  n1 <- ldlasso.obj@n.case
  y <- ldlasso.obj@logOR.norm
 
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
  result <- list(  qp = qp, A = A, r2.vec = r2.vec, b0 = b0 )

  ldlasso.obj@beta <- result$qp$solution[1:p]

  stopifnot(validObject(ldlasso.obj))
  
  return(ldlasso.obj)
  
  }
