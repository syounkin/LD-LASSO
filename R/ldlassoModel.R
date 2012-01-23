ldlassoModel <- function( ldlasso.obj, B = 1e1, alpha = 0.05, tol = 1e-2 ){
  if( is.null(ldlasso.obj@s1) ) {
    ldlasso.obj@s1 <- findS1( ldlasso.obj, B = B, alpha = alpha, tol = tol )
  }
  ldlasso.obj@beta <- solve(ldlasso.obj) 
  return( ldlasso.obj)
}
