printHeaderIdeal <- function(x){
  cat(paste("ideal was called as follows:\n"))
  print(x$call)
  cat("\n")

  cat(paste("Number of Legislators:\t",x$n,"\n"))
  cat(paste("Number of Votes:\t",x$m,"\n"))
  cat(paste("Number of Dimensions:\t",x$d,"\n"))
  cat(paste("Number of Iterations:\t",
            eval(x$call$maxiter,envir=.GlobalEnv),
            "\n"))
  cat(paste("\tThinned By:\t",
            eval(x$call$thin,envir=.GlobalEnv),"\n"))
  cat(paste("\tBurn-in:\t",
            eval(x$call$burnin,envir=.GlobalEnv),"\n\n"))

  invisible(NULL)
}

## summary and print functions
print.ideal <- function(x, ...) {
  if(!inherits(x, "ideal"))
    stop("object passed to print.ideal is not of class ideal\n")

  cat("Markov chain Monte Carlo Analysis of Roll Call Data\n")
  if(x$d==1)
    cat("       (2-parameter item-response modeling)        \n")
  else
    cat("     (multidimensional item-response modeling)     \n")
  cat("===================================================\n\n")

  printHeader(eval(x$call$object))
  printHeaderIdeal(x)

  if(is.null(x$call$file)){
    cat("Ideal Points: Posterior Means\n")
    print(round(x$xbar,2))
    cat("\n")
  }
  invisible(NULL)
}

summary.ideal <- function(object,
                          prob=.95,
                          burnin=NULL,
                          sort=TRUE,
                          include.beta=FALSE,
                          ...){

  if(!inherits(object, "ideal"))
    stop("summary.ideal only defined for objects of class ideal")

  if(!is.null(object$call$file)){
    print(object)
    cat("\n")
    cat(paste("MCMC output was directed to file:",
              object$call$file,
              "\n"))
    cat(paste("no output to summarize in the ideal object",
              match.call()$object,
              "\n"))
    return(invisible(NULL))
  }

  if(is.null(burnin)){
    keep <- checkBurnIn(object,eval(object$call$burnin,envir=.GlobalEnv))
  } else {
    keep <- checkBurnIn(object,burnin)
  }
  
  xm <- NULL
  xsd <- NULL
  bm <- NULL
  bsd <- NULL
  xHDR <- NULL
  bHDR <- NULL
  xResults <- list()
  bResults <- list()
  bSig <- list()

  myHPD <- function(x,prob){
    tmp <- coda::as.mcmc(x)
    return(coda::HPDinterval(tmp,prob))
  }
  
  ## get HPD of x
  xKeep <- object$x[keep,,,drop=FALSE]
  xHDR <- apply(xKeep,
                c(2,3),
                myHPD,
                prob=prob)
  xm <- apply(xKeep,c(2,3),mean)   ## means
  xsd <- apply(xKeep,c(2,3),sd)    ## standard deviations
  dimnames(xHDR)[[1]] <- c("lower","upper") 
  if(length(dim(xHDR))>2)
    xHDR <- aperm(xHDR,c(2,1,3))
  if (length(dim(xHDR))==2)
    xHDR <- t(xHDR)

  ## ################################################################
  ## get beta summaries
  if ((!is.null(object$beta)) && (include.beta)){
    bKeep <- object$beta[keep,,,drop=FALSE]
    bHDR <- apply(bKeep,
                  c(2,3),
                  myHPD,
                  prob=prob)
    dimnames(bHDR)[[1]] <- c("lower","upper")
    if (length(dim(bHDR))>2)
      bHDR <- aperm(bHDR,c(2,1,3))
    
    bm <- apply(bKeep,c(2,3),mean)
    bsd <- apply(bKeep,c(2,3),sd)
    
    ## "significance tests" for discrimination parameters
    ## we have HDR interval of content prob
    sigFunc <- function(x){
      out <- sign(x[1])==sign(x[2])
      labs <- rep(paste(round(prob*100),"% CI",sep=""),2)
      labs[1] <- paste(labs[1],"does NOT overlap 0")
      labs[2] <- paste(labs[2],"overlaps 0")
      out <- factor(out,
                    levels=c("TRUE","FALSE"),
                    labels=labs)
      out
    }
    bSig <- NULL
    bSig <- apply(bHDR,c(1,3),sigFunc)
    bSig <- bSig[,-grep("Difficulty",dimnames(bSig)[[2]])]
  }

  #####################################################################
  ## summarize by party
  pall.final <- NULL
  if(!is.null(object$call$dropList)){
    party <- dropRollCall(eval(object$call$object),
                          eval(object$call$dropList))$legis.data$partyName
    if(is.null(party))
      party <- dropRollCall(eval(object$call$object),
                            eval(object$call$dropList))$legis.data$party
  } else {
    party <- eval(object$call$object)$legis.data$partyName
    if(is.null(party))
      party <- eval(object$call$object)$legis.data$party
  }
  if(!is.null(party)){                       ## we have some party info
    nms <- NULL
    for (b in 1:object$d){       ## loop over dimensions
      pall <- NULL
      pm <- tapply(xm[,b],party,mean)
      pq <- tapply(xm[,b],party,
                   quantile,
                   probs=c(0,prob) + (1-prob)/2)
      for(j in 1:length(pq)){
        pall <- rbind(pall,pq[[j]])
      }
      pall <- cbind(pm,pall)
      pall.final <- rbind(pall.final,pall)
      nms <- c(nms,paste(rownames(pall),": Dimension ",b,sep=""))
    }
    
    colnames(pall.final)[1] <- "Mean"
    rownames(pall.final) <- nms
  }
  
  ## ###################################################################
  ## gather for output
  out <- list(object=match.call()$object,
              xm=xm,xsd=xsd,xHDR=xHDR,
              bm=bm,bsd=bsd,bHDR=bHDR,
              bSig=bSig,
              party.quant=pall.final,
              sort=sort,
              prob=prob)
  
  class(out) <- "summary.ideal"
  
  return(out)
}

print.summary.ideal <- function(x, digits=3, ...){ 
  if (!("summary.ideal" %in% class(x)))
    stop("object passed to print.summary.ideal must be of class summary.ideal")

  cat("Markov chain Monte Carlo Analysis of Roll Call Data\n")
  m <- eval(x$object)$m
  d <- eval(x$object)$d
  if(d==1)
    cat("       (2-parameter item-response modeling)        \n")
  else
    cat("     (multidimensional item-response modeling)     \n")

  printHeader(eval(eval(x$object)$call$object))
  printHeaderIdeal(eval(x$object))
  
  if(!is.null(x$party.quant)) {
    cat("Ideal Points (Posterior Means), by Party\n")
    print(round(x$party.quant,digits))
    cat("\n")
  }  

  ## loop over dimensions
  xResults <- list()
  for(j in 1:d){
    xResults[[j]] <- cbind(x$xm[,j],
                           x$xsd[,j],
                           x$xHDR[,,j])
    cNames <- c("Mean","Std.Dev.","lower","upper")
    dimnames(xResults[[j]])[[2]] <- cNames
    if(x$sort)
      xResults[[j]] <- xResults[[j]][order(xResults[[j]][,1]),]
  }

  for(j in 1:d){
    if(x$sort)
      cat(paste("Ideal Points, Dimension ",j," ",
                "(sorted by posterior means):\n",sep=""))
    else
      cat(paste("Ideal Points, Dimension ",j,":\n",sep=""))
    print(round(xResults[[j]],digits))
    cat("\n")
  }

  ## loop over dimensions
  bResults <- list()
  if(!is.null(x$bm)){
    for(b in 1:(d+1)){
      bResults[[b]] <- cbind(x$bm[,b],
                             x$bsd[,b],
                             x$bHDR[,,b])
      dimnames(bResults[[b]])[[2]][1] <- "Mean"
      dimnames(bResults[[b]])[[2]][2] <- "sd"

      tmpRollCall <- computeMargins(object=eval(eval(x$object)$call$object),
                                    dropList=eval(eval(x$object)$call$dropList))
      
      ## if available, tack on margins data
      if(!is.null(tmpRollCall$voteMargins)){
        bResults[[b]] <- cbind(bResults[[b]],
                               tmpRollCall$voteMargins)
      }
    }
    names(bResults) <- c(paste("Discrimination D",
                               1:d),
                         "Difficulty")
    
    
    ## report statistical tests of significance
    if(length(x$bSig)!=0){
      cat("Statistical tests of discrimination parameters:\n")
      if(d==2){ ## do a cross-tabulation
        cat("dimension 1 (rows) against dimension 2 (columns)\n")
        print(table(x$bSig[[1]],x$bSig[[2]]))
      } else{
        for (j in 1:d){
          cat("Dimension:",j)
          print(table(x$bSig[[j]]))
          cat("\n")
        }
      }
    }
    
    for(j in 1:d){
      cat(paste(names(bResults[[j]]),":\n"))
      theseResults <- bResults[[j]]
      foo <- x$bSig[[j]] == (levels(x$bSig[[j]])[2])
      fooChar <- rep("  ",m)
      fooChar[foo] <- "NS"
      dimnames(theseResults)[[1]] <- paste(dimnames(theseResults)[[1]],
                                           fooChar)
      print(round(theseResults,digits))
      cat("\n")
    }
    cat(paste(names(bResults)[d+1],":\n"))
    print(round(bResults[[d+1]],digits))
  }
  
  cat("\n")
  invisible(NULL)
}


