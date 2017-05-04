## ## for creating a package
.package.Name <- "pscl"

##.First.lib <- function(lib,pkg){
##  library.dynam(.package.Name,
##                pkg,
##                lib)
##}

##.Last.lib <- function(libpath){
##  library.dynam.unload(chname="pscl",libpath=libpath)
##}

.onAttach <- function(...){
  packageStartupMessage("Classes and Methods for R developed in the\n")
  packageStartupMessage("Political Science Computational Laboratory\n")
  packageStartupMessage("Department of Political Science\n")
  packageStartupMessage("Stanford University\n")
  packageStartupMessage("Simon Jackman\n")
  packageStartupMessage("hurdle and zeroinfl functions by Achim Zeileis\n") 
}

.onUnload <- function(libpath){
  library.dynam.unload("pscl",libpath=libpath)
}
