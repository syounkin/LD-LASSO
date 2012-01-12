ldlassoModel <- function( ldlasso.obj ){
  if( is.null(ldlasso.obj@s1) ) {
    ldlasso.obj@s1 <- findS1( ldlasso.obj )
  }
  ldlasso.obj@beta <- solve(ldlasso.obj) 
  return( ldlasso.obj)
}
