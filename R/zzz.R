THISPKG <- "ldlasso"
.ldlassoEnv <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
	version <- packageDescription("ldlasso", field="Version")
	packageStartupMessage(paste("Welcome to LD LASSO version ", version, sep = "" ) )
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
