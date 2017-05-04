## ideal plotting functions
plot.ideal <- function(x, conf.int=0.95, burnin=NULL, ...) {
  if(class(x)!="ideal")
    stop("plot.ideal only available for objects of class ideal")
  
  if(is.null(burnin))
    keep <- checkBurnIn(x,eval(x$call$burnin,envir=.GlobalEnv))
  else
    keep <- checkBurnIn(x,burnin)

  if(x$d == 1)
    plot1d(x,conf.int=conf.int, burnin=burnin, ...)
  else if(x$d==2){
    plot2d(x,burnin=burnin, ...)
  }
  else{ ## more than 2 dimensions    
    xm <- apply(x$x[keep,-1],2,mean,na.rm=T)
    dim <- matrix(rep(NA,(x$n*x$d)),ncol=x$d)
    lbls <- NULL
    for (i in 1:x$d){
      dim[,i] <- xm[seq(i,length(xm),by=x$d)]
      lbls <- c(lbls,paste("Dimension ",i,sep=""))
    }
    pairs(dim,labels=lbls,main="Posterior Mean: Ideal Points")
  }
}

plot1d <- function(x,
                   d=1,
                   conf.int=0.95,
                   burnin=NULL,
                   showAllNames=FALSE,
                   ...){
  if(class(x)!="ideal")
    stop("plot.ideal only available for objects of class ideal")
  
  if(is.null(burnin))
    keep <- checkBurnIn(x,eval(x$call$burnin,envir=.GlobalEnv))
  else
    keep <- checkBurnIn(x,burnin)
  
  checkD(x,d)                    ## check that selected dimension is ok
  checkCI(conf.int)              ## check that confidence interval is ok
  q <- c((1-conf.int)/2, 1-((1-conf.int)/2))  ## quantiles from CI
  
  xm <- x$xbar                   ## xbar
  indx <- order(xm)              ## sort index
  exispar <- par(no.readonly=T)

  myHPD <- function(x,prob){
    tmp <- coda::as.mcmc(x)
    return(coda::HPDinterval(tmp,prob))
  }
  xq <- t(apply(x$x[keep,,1],2,myHPD,prob=conf.int))  ## get HPDs

  ## names etc
  cat(paste("Looking up legislator names and party affiliations\n"))
  cat(paste("in rollcall object",x$call$object,"\n"))
  tmpObject <- dropRollCall(eval(x$call$object),
                            eval(x$call$dropList))
  party <- tmpObject$legis.data$party  ## extract party info
  legis.name <- unclass(dimnames(x$xbar)[[1]])
  longName <- max(nchar(legis.name))
  rm(tmpObject)
  textLoc <- 1.05*min(xq)   ## where to put x labels
  
  if(showAllNames){
    par(mar=c(3,longName*.55,4,2)+0.1,
        oma=rep(0,4))
  } else {
    par(mar=c(3,longName*.75,4,2)+0.1,
        oma=rep(0,4))
  }
  ## title string info
  mainString <- paste("Ideal Points: ",
                      "Posterior Means and ",
                      conf.int*100,
                      "% CIs",sep="")
  if(!is.null(eval(x$call$object)$desc))
    mainString <- paste(eval(x$call$object)$desc,"\n",mainString)
  
  plot(y=c(1-.5,x$n+.5),
       yaxs="i",
       x=1.02*range(xq),
       xaxs="i",
       xlab="",ylab="",
       axes=FALSE,
       type="n",
       ...)
  mtext(mainString,side=3,line=3)
  
  xLims <- 1.01*range(xq)
  
  for (i in 1:x$n){
    if(!showAllNames){
      if((x$n <= 30)||(i %in% as.integer(seq(1,x$n,length=30)))){
        text(x=textLoc,
             y=i,
             labels=legis.name[indx[i]],
             adj=1,xpd=NA)
        lines(x=xLims,y=rep(i,2),lty=1,lwd=.5,col=gray(.45))
      }
    }
    else{
      text(x=textLoc,
           y=i,
           cex=.55,
           labels=legis.name[indx[i]],
           adj=1,xpd=TRUE)
      lines(x=xLims,y=rep(i,2),lty=1,lwd=.5,col=gray(.45))
    }
    
    lines(y=c(i,i),x=xq[indx[i],],lwd=2)
    if (is.null(party)){
      points(y=i,x=xm[indx[i]],col="red",pch=19,xpd=NULL)
    }
    else{
      tbl <- table(party, exclude=NULL)
      cl <- rainbow(length(tbl))
      pt <- xm[indx[i]]
      grp <- match(party[indx[i]],names(tbl))
      points(y=i,
             x=pt,
             pch=19,col=cl[grp],xpd=NULL)    
    }
  }  
  ##par(ps=8)
  ##par(ps=10)
  axis(1)
  axis(3)
  
  par(exispar)
}


plot2d <- function(x,
                   d1=1,
                   d2=2,
                   burnin=NULL,
                   overlayCuttingPlanes=FALSE,
                   ...){
  
  if(class(x)!="ideal")
    stop("plot.ideal only available for objects of class ideal")
  
  if(is.null(burnin))
    keep <- checkBurnIn(x,eval(x$call$burnin,envir=.GlobalEnv))
  else
    keep <- checkBurnIn(x,burnin)

  oCP <- overlayCuttingPlanes   ## local copy of overlayCuttingPlanes
  if(overlayCuttingPlanes){
    if(is.null(x$beta)){
      cat("Item parameters were not stored in ideal object\n")
      cat("No cutting planes will be plotted\n")
      oCP <- FALSE
    }
    if(x$d>2){
      cat("overlay of cutting planes only defined for 2d fits\n")
      oCP <- FALSE
    }
  }
  
  mainString <- paste("Ideal Points:",
                      "Posterior Means")
  if(!is.null(eval(x$call$object)$desc))
    mainString <- paste(eval(x$call$object)$desc,"\n",mainString)
  
  checkD(x,d1)
  checkD(x,d2)
  
  if(d1==d2)
    stop("can't do 2 dimensional summaries of the same dimension\n")
  
  if(is.null(burnin)){    ## use x bar in ideal object
    xm1 <- x$xbar[,d1]
    xm2 <- x$xbar[,d2]
  }
  else{
    xm1 <- apply(x$x[keep,,d1],2,mean)  ## posterior means
    xm2 <- apply(x$x[keep,,d2],2,mean)
  }
  
  if(oCP){
    if(is.null(burnin)){    ## use betabar in ideal object
      b1Bar <- x$betabar[,d1]
      b2Bar <- x$betabar[,d2]
      alphaBar <- x$betabar[,(x$d+1)]
    }
    else{
      bKeep <- x$beta[keep,,,drop=FALSE]
      betaBar <- apply(bKeep,c(2,3),mean)
      b1Bar <- betaBar[,d1]
      b2Bar <- betaBar[,d2]
      alphaBar <- betaBar[,x$d]
    }
  }

  ## get a copy of the rollcall object used by ideal object
  tmpObject <- dropRollCall(eval(x$call$object),
                            eval(x$call$dropList))
  party <- tmpObject$legis.data$party  ## extract party info
  rm(tmpObject)
  if(is.null(party)){
    plot(x=xm1,y=xm2,
         main=mainString,
         type="p",
         xlab=paste("Dimension ",as.character(d1),sep=""),
         ylab=paste("Dimension ",as.character(d2),sep=""),
         xpd=NULL,
         ...)
    if(oCP){
      for(j in 1:x$m)
        abline(a=-alphaBar[j]/b2Bar[j],
               b=-b1Bar[j]/b2Bar[j],
               col=gray(.45))
    }
  }
  else{   ## we have party info
    plot(x=xm1,
         y=xm2,
         main=mainString,
         type="n",
         xlab=paste("Dimension ",as.character(d1),sep=""),
         ylab=paste("Dimension ",as.character(d2),sep=""),
         ...)
    if(overlayCuttingPlanes){
      for(j in 1:x$m)
        abline(a=-alphaBar[j]/b2Bar[j],
               b=-b1Bar[j]/b2Bar[j],
               col=gray(.45),lwd=.5)
    }
    
    tbl <- table(party, exclude=NULL)
    cl <- rainbow(length(tbl))
    for (i in 1:length(tbl)){
      thisParty <- party==names(tbl)[i]
      points(y=xm2[thisParty],
             x=xm1[thisParty],
             pch=16,col=cl[i],
             xpd=NULL)
    }
  }
  
  invisible(NULL)
}

tracex <- function(object,
                   legis=NULL,
                   d=1,
                   conf.int=0.95,
                   multi=FALSE,
                   burnin=NULL,
                   span=.25,
                   legendLoc="topright"){

  warnOption <- options()$warn
  options(warn=-1)
  if(class(object)!="ideal")
    stop("object passed to function tracex must be of class ideal")

  if(is.null(d))
    stop("default value must be supplied for dimension to trace\n")

  if(length(d)>2)
    stop("tracex only works with up to 2 dimensions\n")

  if(!is.character(legis))
    stop("legis must be character (names of legislators)\n")

  Rv <- as.numeric(version$major) + .1*as.numeric(version$minor)
  if(Rv>=2.4)
    old.par <- par(no.readonly=TRUE)  
  else
    old.par <- par()


  ## try matching names
  legis.names <- as.vector(dimnames(object$x)[[2]])
  nLegis <- length(legis)
  p <- list()
  for(j in 1:nLegis){
    p[[j]] <- grep(pattern=paste("^",legis[j],sep=""),
                   x=legis.names)
    if(length(p[[j]])==0){
      cat(paste("could not find legislator",legis[j],"\n"))
      p[[j]] <- NA
    }
    if(!is.null(p[[j]]) & length(p[[j]])>object$d){
      cat(paste("no unique match for legislator",legis[j],"\n"))
      cat("try providing more of the unique identifier for the legislator\n")
      p[[j]] <- NA
    }
  }

  ## process list of matches
  plotName <- rep(NA,length(p))
  for(j in 1:nLegis){    
    if(length(d)==1){
      if(!is.na(p[[j]])){
        p[[j]] <- p[[j]][1]
      }
    }
    if(!is.na(p[[j]])){
      foo <- pmatch(x=legis[j],
                    table=as.vector(dimnames(object$xbar)[[1]]))
      if(!is.na(foo))
        plotName[j] <- dimnames(object$xbar)[[1]][foo]
      cat(paste("matching",
                legis[j],
                "with",
                plotName[j],"\n"))
    }
  }
  p <- p[!is.na(p)]
  plotName <- plotName[!is.na(plotName)]
  names(p) <- plotName
  nLegis <- length(p)

  if(is.null(burnin))
    keep <- checkBurnIn(object,eval(object$call$burnin,envir=.GlobalEnv))
  else
    keep <- checkBurnIn(object,burnin)
  start <- as.numeric(dimnames(object$x)[[1]])[keep][1]
  
  ## #######################################################
  ## one-dimensional stuff
  ## #######################################################
  if(length(d)==1){
    options(warn=0)
    checkD(object,d)

    if(span<=0 | span>=1)
      stop("span must be between 0 and 1")
    
    if((conf.int<=0)||(conf.int>=1))
      stop("conf.int must be between 0 and 1")
    
    rw <- 3
    if (nLegis < 3)
      rw <- nLegis
    count <- 0
    if(rw < 4 & dev.interactive())
      par(mfrow=c(rw,1))
    if(length(legendLoc)==1)
      legendLoc <- rep(legendLoc,nLegis)

    
    for (i in 1:nLegis){
      meat <- object$x[keep,p[[i]],d]
      iter <- as.numeric(dimnames(object$x)[[1]])[keep]
      par(mar=c(4, 4, 4, 2) + 0.1)      

      mainText <- plotName[i]
      if(object$d>1)
        mainText <- paste(mainText,
                          ", Dimension ",d,sep="")
      plot(y=meat,
           x=iter,
           las=1,
           type="l",
           xlab="Iteration",ylab="",
           main=mainText)
      runmean <- cumsum(meat)/1:length(iter)
      lines(iter,runmean,col="red",lwd=3)
      
      lf <- loess(meat~iter,span=span)   ## loess overlay
      lines(iter,predict(lf),col="blue",lwd=3)
      
      xbar <- mean(meat)
      q <- c((1-conf.int)/2, 1-((1-conf.int)/2))
      q <- quantile(meat,q)
      
      abline(h=xbar,lwd=3,col="grey")    ## posterior mean
      abline(h=q[1],lty=2,col="grey")    ## confidence intervals
      abline(h=q[2],lty=2,col="grey")

      count <- count + 1

      ## do legend
      if(!is.null(legendLoc[i]))
        legend(x=legendLoc[i],
               bg="white",
               ncol=1,
               legend=c("Trace",
                 "Cumulative Mean",
                 paste("Moving Average (loess, span=",
                       round(span,2),
                       ")",
                       sep=""),
                 "Posterior Mean",
                 paste(round(100*conf.int),"% Confidence Interval",sep="")),
               lty=c(1,1,1,1,2),
               lwd=c(2,3,3,3,2),
               col=c("black","red","blue","grey","grey"),
               yjust=0,cex=.65)
      
      ## prompt user for more plots if we are in interactive mode
      ##cat("dev.interactive returns ",
      ##    dev.interactive(),"\n")
      ##cat(paste("count=",count,"\n"))
      ##cat(paste("nLegis=",nLegis,"\n"))
      if((count==3) & (nLegis > 3) & (dev.interactive())){
        count <- 0
        readline("Press return/enter to see next set of plots: ")
      }
    }
  }
  
  ## ###################################################################
  ## two-dimensional traceplots
  if(length(d)==2){
    goodD <- d %in% (1:object$d)
    if(!all(goodD))
      stop("invalid dimensions requested in tracex")
    
    col <- rainbow(nLegis)       ## colors
    meat <- list()               ## container for iters to plot
    for(i in 1:nLegis){
      xTraces <- object$x[keep,p[[i]],d[1]]
      yTraces <- object$x[keep,p[[i]],d[2]]
      meat[[i]] <- list(x=xTraces,
                        y=yTraces,
                        col=col[i])
    }
    
    if(!multi){                         ## plot all 2d traces at once
      xRange <- range(unlist(lapply(meat,function(x)x$x)),na.rm=TRUE)
      yRange <- range(unlist(lapply(meat,function(x)x$y)),na.rm=TRUE)

      layout(mat=matrix(c(1,2),1,2,byrow=TRUE),
             widths=c(.7,.3))
      
      par(mar=c(4,4,1,1))
      plot(x=xRange,y=yRange,
           type="n",
           axes=FALSE,
           xlab=paste("Dimension",d[1]),
           ylab=paste("Dimensions",d[2]))
      axis(1,las=1)
      axis(2,las=1)
      lineFunc <- function(obj){
        lines(obj$x,obj$y,col=obj$col)
        points(obj$x[1],obj$y[1],pch=1,col="black",cex=2)
        npoints <- length(obj$x)
        points(obj$x[npoints],obj$y[npoints],
               pch=16,col="black",cex=2)
      }
      
      lapply(meat,lineFunc)
      
      mtext(side=3,outer=FALSE,line=-.5,cex=.75,
            paste("Two-dimensional trace plots, MCMC iterations,\n",
                  eval(object$call$object)$desc,
                  ", Iterations ",
                  start," to ",
                  object$call$maxiter," thinned by ",
                  object$call$thin,sep=""))
      
      ## legend plot, 2nd panel
      par(mar=c(3,0,1,0))
      plot(x=c(0,1),
           y=c(.5,nLegis+.5),
           las=1,
           xlab="",ylab="",xaxs="i",yaxs="i",
           axes=FALSE,type="n")
      ## loop to show lines and legislator name
      for(i in 1:nLegis){
        lines(x=c(0,.15),
              y=rep(i,2),
              lwd=2,
              col=col[i])
        text(x=.25,
             y=i,
             cex=.75,
             plotName[i],
             adj=0)
      }
    }
    
    if(multi){                ## multiple panels, one per legislator
      par(mfrow=c(2,2))
      count <- 0
      for(i in 1:nLegis){
        plot(x=meat[[i]]$x,
             y=meat[[i]]$y,
             type="l",
             las=1,
             xlab=paste("Ideal Point, Dimension ",d[1],sep=""),
             ylab=paste("Ideal Point, Dimension ",d[2],sep=""))
        title(plotName[i])
        count <- count + 1
      }
      if(count==3 & dev.interactive()){
        count <- 0
        readline("Press any key to see next set of plots:  ")
      }
    }
    
  }            ## end 2 dimensional stuff
  
  par(old.par)
  options(warn=warnOption)
  
  invisible(NULL)
}
  
