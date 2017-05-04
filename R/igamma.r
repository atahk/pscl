densigamma <- function(x,alpha,beta){
  if(alpha > 0 & beta > 0 & all(x>0))
    (beta^alpha)/gamma(alpha) * x^(-alpha-1) * exp(-beta/x)
  else
    stop("densigamma: invalid parameters\n")
}

pigamma <- function(q,alpha,beta){
  if(alpha > 0 & beta > 0 & all(q>0))
    1-pgamma(1/q,alpha,beta)
  else
    stop("pigamma: invalid parameters\n")
}

qigamma <- function(p,alpha,beta){
  if(alpha > 0 & beta > 0 & all(p>0) & all(p<1)){
    if((1-p)<=.Machine$double.eps){
      out <- Inf
    }
    else{
      out <- 1/qgamma(1-p,alpha,beta)
    }
  }
  else
    stop("qigamma: invalid parameters\n")
  return(out)
}

rigamma <- function(n,alpha,beta){
  if(alpha > 0 & beta > 0)
    1/rgamma(n=n,alpha,beta)
  else
    stop("rigamma: invalid parameters\n")
}

igammaHDR <- function(alpha,beta,content=.95,debug=FALSE){
  ok <- alpha>0 & beta>0 & content>0 & content<1
  if(!ok)
    stop("igammaHDR: invalid parameters\n")

  func <- function(x0,alpha,beta,content){
    y0 <- densigamma(x0,alpha,beta)
    p0 <- pigamma(x0,alpha,beta)
    p1 <- p0 + content - .Machine$double.eps
    if(p1<1){
      x1 <- qigamma(p0+content,alpha,beta)
      y1 <- densigamma(x1,alpha,beta)
      out <- y0-y1
    }
    else{
       if(debug)
         cat(paste("igammaHDR: upper bound too large",p1,"\n"))
       out <- NaN
     }
    out
  }
 
  tryContent <- content
  flag <- FALSE
  while(!flag){
    if(debug)
      cat(paste("igammaHPR: checking search bounds with content=",
                tryContent,
                "\n"))
    try <- rep(NA,2)
    eps <- .Machine$double.eps
    bounds <- c(eps,qigamma(1-tryContent-eps,alpha,beta))
    try[1] <- func(bounds[1],alpha,beta,content=content)
    try[2] <- func(bounds[2],alpha,beta,content=content)
    if(any(is.nan(try)))
      stop("igammaHPR failed with NaN in func\n")
    if(sign(try[1])!=sign(try[2]))
      flag <- TRUE
    else{
      if(debug){
        cat("igammaHPR: bad bounds\n")
        print(try)
        cat("\n")
      }
      tryContent <- tryContent + .01*(1-tryContent)
    }
  }

  if(debug)
    cat("igammaHPD: done refining search bounds...now optimizing")
  foo <- uniroot(f=func,
                 interval=bounds,
                 tol=1e-12,
                 alpha=alpha,
                 beta=beta,
                 content=content)$root
  if(debug)
    cat("...done\n")

  hpd <- c(foo,
           qigamma(pigamma(foo,alpha,beta)+content,
                   alpha,beta)
           )
  hpd
}

