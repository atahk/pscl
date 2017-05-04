## ideal helper functions

## check validity of a burnin number
## return logical of valid iters
checkBurnIn <- function(object, burnin) {
  theIters <- as.numeric(dimnames(object$x)[[1]])
  if (as.numeric(burnin)>max(theIters))
    stop("burnin greater than number of iterations")
  return (theIters > burnin)
}

checkD <- function(x,d) {
  if ((d<1)||(d>x$d))
    stop("d must be equal to one of the dimensions in the roll call object")
}

checkCI <- function(conf.int) {
  if((conf.int<=0)||(conf.int>=1))
      stop("conf.int must be between 0 and 1")
}

getMean <- function(keep,x){
  xbar <- apply(x[keep,,,drop=FALSE],
                c(2,3),
                mean)
  dimnames(xbar) <- list(dimnames(x)[[2]],
                         dimnames(x)[[3]])
  return(xbar)
}
