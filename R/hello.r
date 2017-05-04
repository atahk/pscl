.onAttach <- function(libname,pkgname){
    ## cat(paste("   pscl",
##               paste(rep(".",floor(getOption("width")*.90 - 4)),collapse=""),
##               "\n",
##               sep="")
##         )

    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2], "\n"))

    #cat("   R classes and methods developed in the\n")
    #cat("   Political Science Computational Laboratory\n")
    #cat("   Department of Political Science, Stanford University\n")
    #cat("   Simon Jackman <jackman@stanford.edu>\n")
    #cat("   http://pscl.stanford.edu\n")
    invisible(NULL)
}

.onUnload <- function(){
    invisible(NULL)
}
