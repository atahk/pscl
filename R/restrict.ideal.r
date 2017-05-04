## function to set restrictions on certain legislators
constrain.legis <- function(obj,
                            dropList=list(codes="notInLegis",lop=0),
                            x, d=1){
  options(warn=-1)
  if(class(obj)!="rollcall")
    stop("object must be of class rollcall")
  
  ## check dimensions of x list items
  if (!is.list(x))
    stop("x must be a list")
  if (length(x)<=d)
    stop("at least d+1 legislators must be constrained")
  options(warn=0)

  ## get working version of roll call object
  cat("constrain.legis: calling dropRollCall to get working version of rollcall object\n")
  rc <- dropRollCall(obj,dropList)
  rc$legis.names <- dimnames(rc$votes)[[1]]
  rc$vote.names <- dimnames(rc$votes)[[2]]
  v <- convertCodes(rc)
  
  n <- nrow(v)
  m <- ncol(v)
  xp <- matrix(rep(0, n*d), nrow=n)
  xpv <- matrix(rep(.01, n*d), nrow=n)
  bp <- matrix(rep(0,m*(d+1)), nrow=m)
  bpv <- matrix(rep(.01, m*(d+1)), nrow=m)

  cat("constrain.legis: generating start values for legislators\n")
  xstart <- x.startvalues(v,d=d)
  xcnst <- xstart*NA
  cat("constrain.legis: implementing constraints\n")
  ## loop over constraints
  for (i in 1:length(x)){
    thisLegis <- names(x)[i]
    if (length(x[[i]])!=d)
      stop("Each element of x must be of length d (dimension of model to be fitted).")
    ind <- pmatch(thisLegis,rc$legis.names)
    if (is.na(ind))
      stop(paste(thisLegis,"was not found in legis.names"))
    cat(paste("matching supplied name",
              thisLegis,"with",
              rc$legis.names[ind],"\n"))
    xp[ind,] <- x[[i]]
    xpv[ind,] <- rep(1e12,d)
    xstart[ind,] <- x[[i]]
    xcnst[ind,] <- x[[i]]
  }
  cat("constrain.legis: re-generating start values for legislators, with constraints\n")
  xstart <- x.startvalues(v,d=d,constraint=xcnst)

  options(warn=-1)
  cat("constrain.legis: generating start values for bill parameters,\n")
  cat("conditional on start values for legislators\n")
  bstart <- b.startvalues(v,xstart,d=d)
  bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                   bstart,
                   bp + 2*sign(bstart-bp)/sqrt(bpv))
  options(warn=0)

  return(list(xp=xp,xpv=xpv,bp=bp,bpv=bpv,x=xstart,b=bstart))
}


constrain.items <- function(obj,
                            dropList=list(codes="notInLegis",lop=0),
                            x, d=1){
  options(warn=-1)
  if(class(obj)!="rollcall")
    stop("object must be of class rollcall")
  options(warn=0)

  rc <- dropRollCall(obj,dropList)
  rc$legis.names <- dimnames(rc$votes)[[1]]
  rc$vote.names <- dimnames(rc$votes)[[2]]
  v <- convertCodes(rc)
  
  n <- nrow(v)
  m <- ncol(v)
  xp <- matrix(rep(0, n*d), nrow=n)
  xpv <- matrix(rep(1, n*d), nrow=n)
  bp <- matrix(rep(0,m*(d+1)), nrow=m)
  bpv <- matrix(rep(0.01, m*(d+1)), nrow=m)

  ## check dimensions of x list items
  if (!is.list(x))
    stop("x must be a list")
  options(warn=-1)
  xstart <- x.startvalues(v,d=d)
  bstart <- b.startvalues(v,xstart,d=d)
  options(warn=0)
  for (i in 1:length(x)) {
    if (length(x[[i]])!=(d))
      stop("Each element of x must be of length d (dimension of model to be fitted).")
    ind <- pmatch(names(x)[i],rc$vote.names)
    if (is.na(ind))
      stop(paste(names(x)[i]," was not found in rc$vote.names"))
     cat(paste("matching supplied name",
              names(x)[i],"with",
              rc$vote.names[ind],"\n"))
    bp[ind,] <- c(x[[i]],0)
    bpv[ind,] <- c(rep(1e12,d),0.01)
    bstart[ind,1:d] <- x[[i]]
  }
  return(list(xp=xp,xpv=xpv,bp=bp,bpv=bpv,xstart=xstart,bstart=bstart))
}

  
