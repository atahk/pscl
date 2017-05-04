betaHPD <- function(alpha,beta,p=.95,plot=FALSE,xlim=NULL,debug=FALSE){
  
  if(is.na(p) | is.nan(p) | p > 1 | p < 0)
    stop("p not between 0 and 1\n")
  
    if(alpha<=1 | beta <=1)
      stop("betaHPD only implemented for alpha and beta both > 1\n")

  ## initialize internal logical flags
  compute <- TRUE
  swap <- FALSE

  if(alpha==beta){
    if(debug)
      cat("symmetric case, alpha=",alpha,"beta=",beta,"\n")
    out <- qbeta((1 + c(-1,1)*p)/2,
                 alpha,beta)
    compute <- FALSE
  }
  if(alpha>beta){
    swap <- TRUE
    alphaStar <- beta
    betaStar <- alpha
  }
  else if(beta>alpha){
    swap <- FALSE
    alphaStar <- alpha
    betaStar <- beta
  }
  if(debug)
    cat("swap=",swap,"\n")
  
  func <- function(x0,alpha,beta){
    y0 <- dbeta(x0,alpha,beta)
    p0 <- pbeta(x0,alpha,beta)
    x1 <- qbeta(p0+p,alpha,beta)
    y1 <- dbeta(x1,alpha,beta)
    out <- abs(y0-y1)
    out
  }
  
  if(compute){
    foo <- try(optimize(f=func,alpha=alphaStar,beta=betaStar,
                        tol=.Machine$double.eps^(.6),
                        interval=c(.Machine$double.eps,
                          qbeta(1-p,
                                alphaStar,betaStar))))
    if(inherits(foo,"try-error")){
      warning("optimization in betaHPD failed\n")
      out <- rep(NA,2)
    }
    else{
      if(debug){
        cat("results of optimization:\n")
        print(foo)
      }
      out <- c(foo$minimum,
               qbeta(pbeta(foo$minimum,alphaStar,betaStar)+p,
                     alphaStar,betaStar)
               )
    }
    if(swap){
      out <- 1-out
      out <- sort(out)
      if(debug){
        cat("swapped back\n")
        print(out)
      }
    }
  }

  ## plotting
  if(plot & all(!is.na(out))){
    xseq <- NULL
    if(length(xlim)==2 & all(!is.na(xlim))){
      if(xlim[2]>xlim[1] & xlim[1] >= 0 & xlim[2] <= 1){
        xseq <- seq(xlim[1]+(.Machine$double.eps^(.25)),
                    xlim[2]+(.Machine$double.eps^(.25)),
                    length=1000)
      }
    }
    if(is.null(xseq))
      xseq <- seq(min(qbeta(.0001,alpha,beta),out[1]),
                  max(qbeta(.9999,alpha,beta),out[2]),
                  length=1000)
    
    plot(xseq,dbeta(xseq,alpha,beta),
         xlab=expression(theta),
         ylab="",
         axes=F,
         type="n")
    axis(1)

    ## get polygon for HDR
    dseq <- seq(out[1],out[2],length=250)
    fx <- dbeta(dseq,alpha,beta)
    polygon(x=c(out[1],dseq,rev(dseq)),
            y=c(0,fx,rep(0,250)),
            border=F,col=gray(.45))
    lines(xseq,dbeta(xseq,alpha,beta))
  }

  out
}
