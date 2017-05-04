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
  packageStartupMessage("Classes and Methods for R developed in the\n",
                        "Political Science Computational Laboratory\n",
                        "Department of Political Science\n",
                        "Stanford University\n",
                        "Simon Jackman\n",
                        "hurdle and zeroinfl functions by Achim Zeileis")
}

.onUnload <- function(libpath){
  library.dynam.unload("pscl",libpath=libpath)
}
