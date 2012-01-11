setMethod("solve", signature( a = "ldlasso", b = "missing" ), function( a, b, ... ) ldlassoSolve(a) )
