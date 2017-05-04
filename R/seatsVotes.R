## seats vote curve class
seatsVotes <- function(x,
                       desc=NULL,
                       method="uniformSwing"){

  xok <- x[!is.na(x)]
  if(length(xok)==0)
    stop("no data to analyze after deleting missings")
  
  if(any(!is.numeric(xok)))
    stop("svCurve only defined for numeric data\n")
  if(any(xok<0))
    stop("negative vote shares not permitted\n")

  if(any(xok>1) & !any(xok>100)){
    cat("proceeding assuming supplied votes are percentages")
    xLocal <- xok/100
  }
  else
    xLocal <- xok

  cl <- match.call()
  if(is.null(cl$method))
    cl$method <- method
  
  m <- 1001
  d0 <- seq(from=-1,to=1,length=m)
  x0 <- rep(NA,m)
  y0 <- rep(NA,m)
  for(i in 1:m){
    x0[i] <- mean(xLocal - d0[i],na.rm=TRUE)
    y0[i] <- mean(xLocal - d0[i] > .5,
                  na.rm=TRUE)
  }

  inBounds <- x0 >= 0 & x0 <= 1
  out <- list(s=y0[inBounds],
              v=x0[inBounds],
              x=xLocal,
              desc=desc,
              call=cl)

  class(out)  <- "seatsVotes"
  out
}

print.seatsVotes <- function(x,...){
  if(!inherits(x,"seatsVotes"))
    cat("print.svCurve only defined for objects of class seatsVotes\n")

  if(is.null(match.call()$digits))
    digits <- .Options$digits
  else
    digits <- match.call()$digits

  if(!is.null(x$desc))
    cat("Seats-Votes Curve:",x$desc,"\n")

  cat("\nSummary of",
      length(x$x),
      "non-missing vote shares:\n")
  print(summary(x$x))

  closestToFifty <- which.min(abs(x$v-.5))
  biasAtFifty <- x$s[closestToFifty] - .5

  cat("Bias at Average Vote Share = .5 is",
      round(biasAtFifty,digits),
      "\n")

  invisible(NULL)
}

summary.seatsVotes <- function(object,...){
  if(!inherits(object,"seatsVotes"))
    cat("summary.svCurve only defined for objects of class seatsVotes\n")
  object
}

plot.seatsVotes <- function(x,
                            type=c("seatsVotes","density"),
                            legend="bottomright",
                            transform=FALSE,
                            ...)
{
  if(!inherits(x,"seatsVotes"))
    cat("plot.svCurve only defined for objects of class seatsVotes\n")

  type <- match.arg(type)
  cl <- match.call()


  ## seats vote curve
  if(type=="seatsVotes"){
    oldpar <- par()
    par(mar=c(4.2,4,5,1),
        las=1)

    if(is.null(cl$xlab))
      xlab <- "Average District Vote"
    else
      xlab <- cl$xlab

    if(is.null(cl$ylab))
      ylab <- "Proportion of Seats Won"
    else
      ylab <- cl$ylab
    
    if(is.null(cl$xlim))
      xlim <- c(0,1)
    else
      xlim <- cl$xlim

    if(is.null(cl$ylim))
      ylim <- c(0,1)
    else
      ylim <- cl$ylim

    if(is.null(cl$xaxs))
      xaxs <- "i"
    else
      xaxs <- cl$xaxs

    if(is.null(cl$yaxs))
      yaxs <- "i"
    else
      yaxs <- cl$yaxs

    plot(x$v,x$s,type="l",
         lwd=3,
         axes=FALSE,
         xaxs=xaxs,
         yaxs=yaxs,
         xlim=xlim,
         ylim=ylim,
         xlab=xlab,
         ylab=ylab,
         ...)
    axis(1,at=seq(0,1,by=.25))
    axis(2,at=seq(0,1,by=.25))
    abline(h=.5,lty=2)
    abline(v=.5,lty=2)
    mtext(side=3,
          at=mean(x$x,na.rm=TRUE),
          "Actual\nResult",cex=.65)
    abline(v=mean(x$x,na.rm=TRUE),lty=2)
    abline(0,1,col=gray(.45),lwd=2)
    
    methodString <- switch(x$call$method,
                           "uniformSwing" = "Uniform Swing")
    
    if(!is.null(x$desc))
      title(paste("Simulated Seats-Votes Curve Using ",methodString,
                  "\n",
                  x$desc,sep=""))
    else
      title(paste("Simulated Seats-Votes Curve Using",methodString))
    
    if(!is.null(legend))
      legend(x=legend,
             col=gray(.45),
             lwd=2,
             lty=1,
             bty="n",
             cex=.65,
             legend="Proportional Representation\n(45 degree line)")
    
    par(oldpar)
  }

  if(type=="density"){
    if(is.null(cl$title)){
      if(is.null(x$desc))
        titleString <- "Density"
      else
        titleString <- paste("Density,",x$desc)
    }
    else
      titleString <- cl$title
    
    if(is.null(cl$xlab))
      xlab <- "Vote Shares"
    else
      xlab <- cl$xlab

    if(is.null(cl$ylab))
      ylab <- ""
    else
      ylab <- cl$ylab

    if(transform){
      transFunc <- function(x){
        v <- log(x/(1-x))
        beta <- sqrt(3)
        xstar <- v*beta
        exp(xstar)/(1+exp(xstar))
      }
      xLocal <- transFunc(x$x)
    }
    else
      xLocal <- x$x
    
    plot(density(xLocal,
                 na.rm=TRUE,
                 from=min(xLocal,na.rm=TRUE),
                 to=max(xLocal,na.rm=TRUE),
                 ...
                 ),
         xlab=xlab,
         ylab=ylab,
         main=titleString,
         axes=FALSE)
    if(transform){
      tcks <- pretty(x$x)
      tcks <- transFunc(tcks)
      axis(1,at=tcks,labels=pretty(x$x))
    }
    else
      axis(1)
    rug(xLocal)
  }
  invisible(NULL)
}
