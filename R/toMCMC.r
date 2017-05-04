## convert ideal object to MCMC object

idealToMCMC <- function(object, burnin=NULL){
  if(class(object)!="ideal")
    stop("idealToMCMC only defined for objects of class ideal")
  
  if(is.null(burnin))
    b <- eval(object$call$burnin)
  keep <- checkBurnIn(object,b)

  iters <- as.numeric(dimnames(object$x[keep,,])[[1]])

  out <- object$x[keep,,]
  if(!is.null(object$beta)){
    J <- dim(object$beta)[3]
    for(j in 1:J){
      thisBeta <- object$beta[keep,,j]
      dimnames(thisBeta)[[2]] <- paste(dimnames(thisBeta[[2]]),
                                       dimnames(object$beta[[3]])[j])
      out <- cbind(out,thisBeta)
    }
  }

  return(coda::mcmc(data=out,
                    start=iters[1],
                    thin=eval(object$call$thin),
                    end=iters[length(iters)]
                    )
         )
}
