#setMethod( "initialize", "ldlasso", function( .Object, ...){
#  callNextMethod( .Object, ... )
#})

#setMethod( "initialize", "ldlassoSolution", function( .Object, ... ){
#  .Object <- callNextMethod()
#  .Object@beta <- ldlassoSolve(.Object)
#  .Object
#})

