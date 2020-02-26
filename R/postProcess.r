## post-process an ideal object
postProcess <- function(object,
                        constraints="normalize",
                        debug=FALSE){
  if(!inherits(object, "ideal"))
    stop("postProcess only defined for objects of class ideal")

  ## not a list of normalizing constraint (i.e., usually "normalize=TRUE" with d=1)
  if(!is.list(constraints)){
    if(constraints=="normalize"){
      ## get constraints needed for mean zero, standard deviation one restriction
      tMat <- getNormalizingTransform(object)  ## coefficients for linear map
      if(debug){
        cat("transformation matrix is:\n")
        print(tMat)
      }
      newObject <- implementConstraints(object,tMat,debug)
    }
  }
  
  if(is.list(constraints))
    newObject <- postProcessAffine(object,constraints,debug)

  return(newObject)
}

getNormalizingTransform <- function(object){
  n <- object$n
  m <- object$m
  d <- object$d

  offsets <- apply(object$xbar,2,mean)  
  s <- apply(object$xbar,2,sd)

  coefs <- 1/s * diag(d)               ## d by d
  coefs <- rbind(coefs,-offsets/s)     ## d+1 by d
  return(coefs)
}

affineTrans <- function(x,target){
  d <- dim(x)[2]
  x0 <- cbind(x,1)
  if(d>1){
    zeroMat <- 0*x0
    A <- rbind(cbind(x0,zeroMat),
               cbind(zeroMat,x0))
    b <- as.vector(target)
  }
  if(d==1){
    A <- x0
    b <- target
  }
  foo <- solve(A)%*%b
  foo <- matrix(foo,nrow=d+1)
  return(foo)
}

postProcessAffine <- function(object,constraints,debug){
  d <- object$d
  n <- object$n
  m <- object$m
  nSavedIters <- dim(object$x)[1]
  theIters <- dimnames(object$x)[[1]]
  keep <- checkBurnIn(object,
                      burnin=object$call$burnin)
  
  nCon <- length(constraints)
  if(nCon != (d+1)){
    cat("postProcess is currently only implements as many constraints\n")
    cat("as there are dimensions plus one.\n")
    stop()
  }
  lengthCon <- lapply(constraints,length)
  if(any(lengthCon != d))
    stop("each constraint must have the same number of dimensions as the fitted model")
  
  ## form target matrix
  target <- matrix(NA,d+1,d)
  for(i in 1:(d+1))
    target[i,] <- constraints[[i]]
  if(debug){
    cat("target:\n")
    print(target)
  }
    
  ## form id vector, where are the named legislators in the ideal object?
  legis.names <- dimnames(eval(object$call$object)$votes)[[1]]
  if(is.null(legis.names)){
    cat("can not find legislator names to match against\n")
    cat(paste("either the original roll call object",
              object$call$object,
              "has been deleted\n"))
    cat("or the vote component of the roll call object has been deleted?\n")
    cat("terminating postProcess with an error\n")
    stop()
  }
  
  legis <- names(constraints)
  ind <- rep(NA,d+1) 
  for(i in 1:nCon){
    p <- grep(pattern=paste("^",legis[i],sep=""),
              x=legis.names)
    if(length(p)==0)
      stop("could not find the named legislator in the rollcall object")
    else
      ind[i] <- p
  }
  cat(paste("matching legislators",legis,"\n"))
  
  ## initialize output objects
  newX <- NA * object$x
  dimnames(newX) <- dimnames(object$x)
  
  haveBeta <- eval(object$call$store.item)
  if(haveBeta){
    cat("will also transform item/bill parameters\n")
    newBeta <- NA * object$beta
    dimnames(newBeta) <- dimnames(object$beta)
  }
  
  ## now loop over iterations
  for(iter in 1:nSavedIters){
    cat(paste("post-processing iteration",theIters[iter],"\n"))

    x0 <- object$x[iter,ind,,drop=TRUE]
    x0 <- matrix(x0,d+1,d)
    tMat <- affineTrans(x0,target=target)
    if(debug){
      cat("transformation matrix:\n")
      print(tMat)
    }
    thisX <- cbind(object$x[iter,,],1)
    tmpX <- thisX%*%tMat
    newX[iter,,] <- tmpX
    
    ## now transform beta (and alpha), if available
    if(haveBeta){
      tMatStar <- rbind(t(tMat),
                        c(rep(0,d),1))
      itMat <- try(solve(tMatStar))
      if(!inherits(itMat,"try-error")){
        itMat[1:d,d+1] <- -itMat[1:d,d+1]   ## sign fix for minus intercept
        beta0 <- object$beta[iter,,]
        tmpBeta <- beta0%*%itMat
        newBeta[iter,,] <- tmpBeta
        if(debug){
          cat("inverse transformation matrix:\n")
          print(itMat)
        }
      }
      
      if(debug){
        muPP <- pnorm(cbind(newX[iter,,],-1)%*%t(tmpBeta))
        mu <- pnorm(cbind(object$x[iter,,],-1)%*%t(beta0))
        cat("sanity check, comparison of predictions from original and post-processed:\n")
        print(summary(as.vector(mu-muPP)))
      }
    } 
  }
  
  ## new ideal object
  newObject <- object

  ## new ideal point samples
  newObject$x <- newX
  dimnames(newObject$x) <- dimnames(object$x)
  ## ideal point posterior means
  newObject$xbar <- getMean(keep,newObject$x)
  
  ## for beta?
  if(haveBeta){
    newObject$beta <- newBeta
    dimnames(newObject$beta) <- dimnames(object$beta)
    newObject$betabar <- getMean(keep,newObject$beta)
  }
  
  return(newObject)
}


## implementConstraints
implementConstraints <- function(object,tMat,debug){
  haveBeta <- eval(object$call$store.item)
  if(haveBeta){
    d <- dim(tMat)[2]
    tMatStar <- rbind(t(tMat),
                      c(rep(0,d),1))
    itMat <- try(solve(tMatStar))
    if(inherits(itMat,"try-error"))
      stop("could not compute normalizing transformation for item parameters\n") 
    newBeta <- array(NA,dim(object$beta))
  }

  ## get burnin
  ##if(is.symbol(object$call$burnin)){
  ##  burnin <- eval(object$call$burnin)
  ##} else {
  ##  burnin <- object$call$burnin
  ##}
  keep <- checkBurnIn(object,
                      burnin=eval(object$call$burnin))
  
  nSavedIters <- dim(object$x)[1]
  newX <- array(NA,dim(object$x))
  newObject <- object ## copy ideal object
  
  ## loop over iterations, implementing transformation
  for(iter in 1:nSavedIters){
    thisX <- cbind(newObject$x[iter,,],1)  ## add intercept for translation
    newX[iter,,] <- thisX%*%tMat           ## transformation
    
    ## now transform beta (and alpha), if available
    if(haveBeta){
      beta0 <- object$beta[iter,,]
      tmpBeta <- beta0%*%itMat
      newBeta[iter,,] <- tmpBeta
      
      if(debug){
        muPP <- pnorm(cbind(newX[iter,,],-1)%*%t(tmpBeta))
        mu <- pnorm(cbind(object$x[iter,,],-1)%*%t(beta0))
        cat("sanity check, comparison of predictions from original and post-processed:\n")
        print(summary(as.vector(mu-muPP)))
      }
      
    }
  }

  ## gather up for new ideal object
  newObject$x <- newX
  dimnames(newObject$x) <- dimnames(object$x)
  ## new posterior means
  newObject$xbar <- getMean(keep,newObject$x)
  
  ## for Beta?
  if(haveBeta){
    newObject$beta <- newBeta
    dimnames(newObject$beta) <- dimnames(object$beta)
    newObject$betabar <- getMean(keep,newObject$beta)
  }

  return(newObject)
}
