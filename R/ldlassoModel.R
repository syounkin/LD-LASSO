ldlassoModel <- function( ldlasso.obj ){
  if( is.null(ldlasso.obj@s1) ) {
    ldlasso.obj@s1 <- findS1( ldlasso.obj )
  }
  return( ldlasso.obj)
}

findS1 <- function( ldlasso.obj ){
  # place holder for function to find s1 given s2
  cat( "Null value for s1 so setting s1 = 1\n\n")
  return(1)
}
