#setValidity("ldlassoSolution", function(object){
#   if( any( beta(object)[indexmat(object)[,1]] < 0 ) ) 
#     return("LD LASSO constraint violated")
#})
