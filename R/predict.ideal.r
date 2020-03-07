## predict method for class ideal
predict.ideal <- function(object,
                          cutoff=0.5,
                          burnin=NULL,
                          ...) {
  if(!inherits(object, "ideal"))
    stop("predict.ideal only defined for objects of class ideal\n")
  
  if(is.null(object$beta)){
    cat("Beta values must have be stored in ideal object to make predictions")
    stop("try re-fitting model with store.item=TRUE.")
  }

  if(is.null(burnin))
    keep <- checkBurnIn(object,eval(object$call$burnin))
  else
    keep <- checkBurnIn(object,burnin)  ## check that start is valid

  ## get votes into shape for prediction from ideal
  cat(paste("predict.ideal: Working with rollcall object",
            object$call$object,
            "\n"))
  rcObj <- try(dropRollCall(eval(object$call$object),
                            eval(object$call$dropList)),
               silent=TRUE)
  if(inherits(rcObj,"try-error") | is.null(rcObj$votes)){
    cat(paste("The ideal object ",
               as.name(object),
              " was fitted using the\n",
              "rollcall object" ,
              as.name(object$call$object),
              " which can no longer be found,\n",
              "or does not have a votes matrix as one of its components.",
              sep=""))
    stop("Prediction can not proceed.")
  }
  
  if(checkVotes(rcObj))  ## codes ok, should be ok
    stop("bad votes in rollcall object, can't generate predictions")
  cat("\n")
  votes <- convertCodes(rcObj)       ## convert to 0, 1 and NAs 
  
  ## predictions at posterior means
  predprob <- matrix(NA, ncol=ncol(votes), nrow=nrow(votes))
  dimnames(predprob) <- dimnames(votes)
  pred <- predprob
  correct <- predprob
  if(!is.null(burnin)){
    cat("Computing posterior means using ideal object.\n")
    x1 <- matrix(apply(object$x[keep,],
                       c(2,3),
                       mean),
                 nrow=object$n,
                 ncol=object$d,
                 byrow=TRUE)
    x1 <- cbind(x1,-1)            ## negative intercept !!! SDJ 05/15/07
    b <- matrix(apply(object$beta[keep,-1],2,mean),
                nrow=object$m,ncol=object$d+1,byrow=TRUE)
  }
  else{
    cat("Using posterior means in ideal object.\n")
    x1 <- cbind(object$xbar,-1.0)   ## negative intercept !!! SDJ 01/22/07
    b <- object$betabar
  }
  mu <- tcrossprod(x1,b)        ## this should be n by (d+1) times (d+1) by m
  predprob <- pnorm(mu)
  pred <- predprob >= cutoff
  correct <- votes==pred
  correct[is.na(votes)] <- NA

  lp <- apply(correct,1,tally)*100   ## legislator-specific
  pp <- NULL                         ## by party
  party <- eval(object$call$object)$legis.data$party
  if(!is.null(party))
    pp <- tapply(lp,party,mean)
  
  out <- list(pred.probs=predprob,
              prediction=pred,
              correct=correct,
              legis.percent=lp,
              vote.percent=apply(correct,2,tally)*100,
              yea.percent=(sum(correct[votes==1],na.rm=T)/
                           sum(!is.na(correct[votes==1])))*100,
              nay.percent=(sum(correct[votes==0],na.rm=T)/
                           sum(!is.na(correct[votes==0])))*100,
              party.percent=pp,
              overall.percent=(sum(correct,na.rm=T)/sum(!is.na(correct)))*100,
              ideal=match.call()$object)
  class(out) <- "predict.ideal"
  out
}

tally <- function(x){
  sum(x,na.rm=T)/sum(!is.na(x))
}

print.predict.ideal <- function(x,digits=2,...) {
  cat(paste("Predictions using ideal object",
            x$ideal,
            "\n"))
  cat(paste(x$ideal,"uses rollcall object",eval(x$ideal)$call$object,"\n"))
  rcObj <- eval(x$ideal)$call$object
  printHeader(eval(rcObj))
            
  cat("Predictions calculated using posterior means for\n")
  cat("legislators ideal points and bill parameters\n\n")
  cat("Percent correctly predicted:\n")
  cat(paste("\tOverall:\t",round(x$overall.percent,digits),"%\n",sep=""))
  cat(paste("\tYeas:\t\t",round(x$yea.percent,digits),"%\n",sep=""))
  cat(paste("\tNays:\t\t",round(x$nay.percent,digits),"%\n\n",sep=""))
  cat("Percent Correctly Predicted by Legislator\n")
  mat <- round(as.matrix(x$legis.percent),digits)
  colnames(mat) <- "Percent"
  print(mat)
  if(!is.null(x$party.percent)) {
    cat("\nPercent Correctly Predicted by Legislator, Party Average\n")
    print(round(x$party.percent,digits))
  }
  mat <- round(as.matrix(x$vote.percent),digits)
  colnames(mat) <- "Percent"
  cat("\nPercent Correctly Predicted by Vote\n")
  print(mat)
  cat("\n")
  invisible(NULL)
}

plot.predict.ideal <- function(x,
                               type=c("legis","votes"),
                               ...){
  if(!inherits(x, "predict.ideal"))
    stop("plot.predict.ideal only defined for objects of class predict.ideal")

  localType <- match.arg(type)

  d <- eval(x$ideal)$d
  desc <- eval(eval(x$ideal)$call$object)$desc
  rcObj <- extractRollCallObject(eval(x$ideal))
  
  if(localType=="legis"){
    ## plot percent correctly predicted against posterior mean of ideal point
    ## dimension by dimension
    xbar <- eval(x$ideal)$xbar
    n <- eval(x$ideal)$n
    
    ## party colors
    party <- rcObj$legis.data$party
    if(!is.null(party)){
      tbl <- table(party,exclude=NULL)
      cl <- rainbow(length(tbl))
      grp <- match(party,names(tbl))
      col <- cl[grp]
    }
    else
      col <- rep("black",n)

    for(b in 1:d){
      plot(y=x$legis.percent,
           x=xbar[,b],
           col=col,
           xlab="Ideal Point (Posterior Mean)",
           ylab="Voting Decisions Correctly Predicted (%)",
           pch=16)
      titleString <- paste(desc,"\n",
                           "Percent Correctly Predicted, by Legislator's Ideal Point")  
      if(d>1)
        titleString <- paste(titleString,
                             "(Dimension",b,")")
      title(titleString)
    }
  }

  #########################################################################
  ## type = votes
  #########################################################################
  if(localType=="votes"){
    if(is.null(rcObj$voteMargins))
      rcObj <- computeMargins(rcObj)
    margin <- rcObj$voteMargins
    
    if(is.null(margin))
      stop("failed to find or computes votes margins")
    
    ## margin information
    vote.percent <- t(apply(margin,1,function(x)x/sum(x[1:2])*100))
    margin <- cbind(margin,vote.percent[,1:2])
    
    #dimnames(margin)[[2]] <- c("Yea","Nay","NA",
    #                           "Yea (proportion of those voting)",
    #                           "Nay (proportion of those voting)")
    plot(x=jitter(margin[,4]),   ## percent voting yes
         y=x$vote.percent,
         xlab="Losing Coalition (%, excluding NAs, jittered)",
         ylab="Voting Decisions Correctly Predicted (%, excluding NA)")
    titleString <- paste(desc,"\n",
                         "Percent Correctly Predicted, by Vote Margin")
    title(titleString)
  }

  invisible(NULL)
}
