## IDEAL
ideal <- function(object,
                  codes=object$codes,
                  dropList=list(codes="notInLegis",lop=0),
                  d=1,
                  maxiter=10000,
                  thin=100,
                  burnin=5000,
                  impute=FALSE,
                  normalize=FALSE,
                  meanzero=normalize,
                  priors=NULL,
                  startvals="eigen",
                  store.item=FALSE,
                  file=NULL,
                  verbose=FALSE,
                  use.voter=NULL){

  cat("ideal: analysis of roll call data via Markov chain Monte Carlo methods.\n\n")

  ## calling args, some evaluated if symbols, for future use
  cl <- match.call()
  
  if(is.null(cl$d) | is.symbol(cl$d))
    cl$d <- eval(d,parent.frame())
  if(is.null(cl$codes))
    cl$codes <- codes
  if(is.null(cl$dropList))
    cl$dropList <- dropList
  if(is.null(cl$maxiter) | is.symbol(cl$maxiter))
    cl$maxiter <- eval(maxiter,parent.frame())
  if(is.null(cl$thin) | is.symbol(cl$thin))
    cl$thin <- eval(thin,parent.frame())
  if(is.null(cl$burnin) | is.symbol(cl$burnin))
    cl$burnin <- eval(burnin,parent.frame())
  if(is.null(cl$impute))
    cl$impute <- impute
  if(is.null(cl$store.item) | is.symbol(cl$store.item))
    cl$store.item <- eval(store.item,parent.frame())
  if(is.null(cl$normalize))
    cl$normalize <- normalize
  if(is.null(cl$verbose))
    cl$verbose <- verbose
  
  mda <- FALSE
  
  ## check validity of user arguments
  if (!("rollcall" %in% class(object)))
    stop("object must be of class rollcall")
  if(((d%%1) != 0) || (d<1)){
    stop("d is not a positive integer")
  }

  if(((thin%%1)!=0) || (thin<1)) {
    stop("thin is not a positive integer")
  }

  if(((maxiter%%1)!=0) || (maxiter<1)) {
    stop("maxiter is not a positive integer")
  }

  if(!is.list(dropList))
    stop("dropList must be a list")

  if(!is.list(codes))
    stop("codes must be a list")

  ##check iterations and thinning
  if ((maxiter%%thin)!=0) {
    stop("maxiter must be a multiple of thin")
  }

  if ((burnin%%thin)!=0) {
    stop("burnin must be a multiple of thin")
  }

  if (burnin >= maxiter)
    stop("burnin must be less than maxiter")

  if(!is.null(normalize) & d>1){
    cat("normalize option is only meaningful when d=1\n")
  }

  if(normalize != meanzero){
    normalize <- meanzero
    cat("meanzero option is being phased out; normalize provides the same functionality\n")
    cat(paste("For now, we will use your supplied value of meanzero, proceeding with normalize=",
              meanzero,"\n"))
  }
              
  ## pre-process rollcall object
  tmpObject <- object
  if(!is.null(codes)){
    tmpObject$codes <- codes
    if(checkCodes(tmpObject$codes))
      stop("supplied codes fail redundancy checks")
  }
  if(!is.null(dropList)){
    if(verbose)
      cat(paste("Subsetting rollcall object",
                as.name(cl$object),
                "using dropList\n"))
    y <- dropRollCall(tmpObject,dropList)  ## any subsetting to do?
  }
  else
    y <- tmpObject
  rm(tmpObject)
  
  n <- dim(y$votes)[1]
  m <- dim(y$votes)[2]
  legis.names <- dimnames(y$votes)[[1]]
  vote.names <- dimnames(y$votes)[[2]]

  if (!is.null(use.voter)) {
      if (!is.vector(use.voter))
          stop("use.voter must be a vector of length n")
      if (n != length(use.voter))
          stop("use.voter must be a vector of length n")
  }

  ## map roll call votes into binary format required by ideal
  if(verbose){
    printCodes(codes)
    cat("\n")
  }
  
  if(checkVotes(y$votes,codes))
    stop("rollcall: can't map all votes using supplied codes")

  v <- convertCodes(y,codes)  ## convert to zeros and ones and NAs
  
  ## using a file for storage
  usefile <- !is.null(file)      

  ## check to see how much information will need to be stored
  numrec <- (maxiter-burnin)/thin+1
  if (interactive() & verbose &
      ((store.item)&&((n+m)*d*numrec>2000000))
      ||
      ((!store.item)&&((n*d*numrec)>2000000))
      ){
    ans <- readline(paste("The current call to ideal will result in a large object that\n",
                          "will take up a large amount of memory.  Do you want to\n",
                          "continue with the current configuation? (y/n): ",
                          sep=""))
    
    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated execution of ideal.")
  }

  if (interactive() & verbose & numrec>1000) {
    ans <- readline(paste("You are attempting to save ",numrec," iterations.  This\n",
                          "could result in a very large object and cause memory problems.\n",
                         "Do you want to continue with the current call to ideal? (y/n): ",
                         sep=""))
    
    if ((substr(ans, 1, 1) == "n")||(substr(ans, 1, 1) == "N"))
      stop("User terminated execution of ideal.")
  }

  cat(paste("Ideal Point Estimation\n\nNumber of Legislators\t\t",
              n,"\nNumber of Items\t\t\t", m, "\n\n"))
  xp <- xpv <- bp <- bpv <- NULL

  ####################################################################
  ## check priors
  ####################################################################
  if(verbose)
    cat("checking for any user-supplied priors...\n")
  if(!is.null(priors)){
    if(!is.list(priors))
      stop("priors must be a list")
    
    if(all(unlist(lapply(priors,is.null))))
      stop("priors supplied in a list, but all elements are NULL")
    
    if(sum(unlist(lapply(priors,is.na)))>0)
      stop("priors contain missing values, which is not allowed")
  
    ## now check individual elements of prior list
    if(!is.null(priors$xp)){
      if(length(priors$xp)==1)     ## user supplied a scalar
        xp <- matrix(priors$xp,n,d)
      ## coerce a vector to a matrix
      if(length(priors$xp)>1 & d==1 & !is.matrix(priors$xp))
        xp <- matrix(priors$xp,n,d)
      if(is.matrix(priors$xp))
        xp <- priors$xp
    }
    else{
      if(verbose)
        cat("no prior means supplied for ideal points,\n",
            "setting to default of 0\n")
      xp <- matrix(0,n,d)
    }
    
    if(!is.null(priors$xpv)){
      if(length(priors$xpv)==1)    ## user supplied a scalar
        xpv <- matrix(priors$xpv,n,d)
      ## coerce a vector to a matrix
      if(length(priors$xpv)>1 & d==1 & !is.matrix(priors$xpv))
        xpv <- matrix(priors$xpv,n,d)
      if(is.matrix(priors$xpv))
        xpv <- priors$xpv
    } else {
      if(verbose)
        cat("no prior precisions supplied for ideal points,\n",
            "setting to default of 1\n")
      xpv <- matrix(1,n,d)
    }
    
    if(!is.null(priors$bp)){
      if(length(priors$bp)==1)    ## user supplied a scalar
        bp <- matrix(priors$bp,m,d+1)
      if(is.matrix(priors$bp))
        bp <- priors$bp
    } else {
      if(verbose)
        cat("no prior means supplied for item parameters,\n",
            "setting to default to 0\n")
      bp <- matrix(0,m,d+1)
    }
    
    if(!is.null(priors$bpv)){
      if(length(priors$bpv)==1){   ## user supplied a scalar
        bpv <- matrix(priors$bpv,m,d+1)
      }
      if(is.matrix(priors$bpv)){
        bpv <- priors$bpv
      }
    } else {
      if(verbose){
        cat("no prior precisions supplied for item parameters,\n",
            "setting to default of .04\n")
      }
      bpv <- matrix(.04,m,d+1)
    }
    
    if (((nrow(xp) != n)||(ncol(xp) != d)) || ((nrow(xpv)!=n)||(ncol(xpv)!=d))) {
      stop("Dimensions of xp or xpv not n by d")
    }
    
    if (((nrow(bp) != m)||(ncol(bp) != (d+1))) || ((nrow(bpv)!=m)||(ncol(bpv)!=(d+1)))) {
      stop("Dimensions of bp or bpv not m by d+1")
    }
  }
  
  ## ##################################################################
  ## if we get this far with priors still NULL
  ## then revert to defaults
  ## ##################################################################
  if(is.null(xp)){
    if(verbose)
      cat("setting prior means for ideal points to all zeros\n")
    xp <- matrix(0,n,d)
  }
  if(is.null(xpv)){
    if(verbose)
      cat("setting prior precisions for ideal points to all 1\n")
    xpv <- matrix(1,n,d)
  }
  if(is.null(bp)){
    if(verbose)
      cat("setting prior means for item parameters to all zeros\n")
    bp <- matrix(0,m,d+1)
  }
  if(is.null(bpv)){
    if(verbose)
      cat("setting prior precisions for item parameters to all 0.04\n")
    bpv <- matrix(0.04,m,d+1)
  }
  
  xp <- as.vector(xp)
  xpv <- as.vector(xpv)
  bp <- as.vector(bp)
  bpv <- as.vector(bpv) 

  ################################################################
  ## check for start values - create if not supplied
  ################################################################
  if(verbose)
    cat("\nchecking start values...\n")
  xstart <- NULL
  bstart <- NULL
  options(warn=-1)

  if(!is.list(startvals)){
    if(startvals=="eigen" | is.null(startvals)){
      xstart <- x.startvalues(v,d=d,verbose=verbose)
      bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
      bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                       bstart,
                       bp + 2*sign(bstart-bp)/sqrt(bpv))
    }

    if(startvals=="random"){
      if(verbose)
        cat("generating start values for ideal points by iid sampling from N(0,1)\n")
      xstart <- matrix(rnorm(n*d),n,d)
      bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
      bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                       bstart,
                       bp + 2*sign(bstart-bp)/sqrt(bpv))
    }
  }

  ## user has passed something in startvals
  if(is.list(startvals)){
    cat("found user-supplied list in startvals\n")
    cat("starvals is a list containing:\n")
    print(names(startvals))
    
    if(!is.null(startvals$x)){
      if(length(startvals$x) != n*d)
        stop("supplied start values for x is not n by d")
      if(d==1)
        xstart <- matrix(startvals$x,ncol=1)
      else
        xstart <- startvals$x
      if (sum(is.na(xstart))!=0)
        stop("xstart contains missing values")
    }

    if(!is.null(startvals$b)){
      if(length(startvals$b) != m*(d+1))
        stop("length of bstart not m by d+1")
      bstart <- startvals$b
      if(sum(is.na(bstart))!=0)
        stop("bstart contains missing values")
    } 
  } 

  ## final check
  if(is.null(xstart)){
    cat("no user-supplied start values found\n")
    xstart <- x.startvalues(v,d,verbose=TRUE)
  }
  if(is.null(bstart)){
    bstart <- b.startvalues(v,xstart,d=d,verbose=verbose)
    bstart <- ifelse(abs(bstart - bp) < 2/sqrt(bpv),
                     bstart,
                     bp + 2*sign(bstart-bp)/sqrt(bpv))
  }

  ## report to user
  if(verbose){
    if(n<501){
      cat("using the following start values for ideal points:\n")      
      print(xstart)
    } else {
      cat("using the following start values for ideal points (summary follows):\n")      
      print(summary(xstart))
    }

    if(m<501){
      cat("using the following start values for item parameters:\n")
      print(bstart)
    } else {
      cat("using the following start values for item parameters (summary follows):\n")
      print(summary(bstart))
    }
  }

  xstart <- as.vector(xstart)
  bstart <- as.vector(bstart)
  
  options(warn=0)

  ##############################################################
  ## end error checking
  ##############################################################

  yToC <- ifelse(is.na(v), 9, v)
  yToC <- as.vector(yToC)
  cat("\nStarting MCMC Iterations...\n")

  ## ############################################
  ## two versions, one with usefile option
  ## ############################################
  if (usefile) {
    if (length(legis.names) == n) {
      cat(paste("\"",c("Iteration",legis.names),"\"", sep="", collapse=","),
          file=file)
    }
    else {
      cat(paste("\"",c("Iteration",paste("x", 1:n, sep="")),"\"",
                sep="", collapse=","),
          file=file)
    }
    
    if (store.item){
      cat(",",
          paste("\"",
                c(paste("b",
                        as.vector(apply(expand.grid(1:m,1:(d+1)),1,paste,collapse=".")),
                        sep=".")),"\"",
                sep="", collapse=","),
          sep="", file=file, append=TRUE)
    }
    cat("\n", file=file, append=TRUE)
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n),           #1
                 as.integer(m),           #2
                 as.integer(d),           #3
                 as.double(yToC),         #4
                 as.integer(maxiter),     #5
                 as.integer(thin),        #6
                 as.integer(impute),      #7
                 as.integer(mda),         #8
                 as.double(xp),           #9
                 as.double(xpv),          #10
                 as.double(bp),           #11
                 as.double(bpv),          #12
                 as.double(xstart),       #13
                 as.double(bstart),       #14
                 xoutput=0,               #15
                 boutput=0,               #16
                 as.integer(burnin),      #17
                 usefile,                 #18
                 as.logical(store.item),              #19
                 as.character(file),      #20
                 as.logical(verbose),     #21
                 as.logical(!is.null(use.voter)),   #22 
                 as.integer(use.voter))   #23
  }
  ## not saving output to file, saving output to memory
  else if (!store.item) {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), 
                 as.integer(m), 
                 as.integer(d), 
                 as.double(yToC), 
                 as.integer(maxiter), 
                 as.integer(thin), 
                 as.integer(impute),
                 as.integer(mda),
                 as.double(xp), 
                 as.double(xpv), 
                 as.double(bp),
                 as.double(bpv), 
                 as.double(xstart), 
                 as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=as.double(rep(0,m*(d+1)*numrec)),
                 as.integer(burnin),
                 usefile, 
                 as.logical(store.item), 
                 as.character(file),
                 as.logical(verbose), 
                 as.logical(!is.null(use.voter)), 
                 as.integer(use.voter))
  }
  else {
    output <- .C("IDEAL",
                 PACKAGE=.package.Name,
                 as.integer(n), as.integer(m), as.integer(d), as.double(yToC),          
                 as.integer(maxiter), as.integer(thin), as.integer(impute),
                 as.integer(mda),
                 as.double(xp), as.double(xpv), as.double(bp),
                 as.double(bpv), as.double(xstart), as.double(bstart),
                 xoutput=as.double(rep(0,n*d*numrec)),
                 boutput=as.double(rep(0,m*(d+1)*numrec)), 
                 as.integer(burnin),
                 usefile, 
                 as.logical(store.item), 
                 as.character(file),
                 as.logical(verbose), 
                 as.logical(!is.null(use.voter)), 
                 as.integer(use.voter))
  }

  cat("\n")

  ## parse returns from C
  xbar <- NULL
  betabar <- NULL
  if (!usefile) {
    itervec <- seq(burnin,maxiter,by=thin)
    keep <- itervec > burnin

    ## ideal points
    x <- array(output$xoutput,
               c(n,d,numrec))
    
    ## reshape to iteration first format
    x <- aperm(x,c(3,1,2))   
    dimnames(x) <- list(itervec,
                        legis.names,
                        paste("D",1:d,sep=""))
    if(verbose)
      cat("...computing posterior means for ideal points...")
    xbar <- getMean(keep,x)
    if(verbose)
      cat("done\n")

    ###############################################################
    ## item parameters
    if(store.item){
      b <- array(output$boutput,c(m,d+1,numrec))  ## votes by parameters by iters
      dimnames(b) <- list(vote.names,
                          c(paste("Discrimination D",1:d,sep=""),
                            "Difficulty"),
                          itervec)
      ## reshape to iteration first format
      b <- aperm(b,c(3,1,2))                      ## iters by votes by parameters
      if(verbose)
      cat("...computing posterior means for item parameters...")
      betabar <- getMean(keep,b)
      if(verbose)
        cat("done\n")
    } else {
      b <- NULL
    }
  } else {   ## output went to a file
    b <- x <- NULL
  }
 
  ## wrap up for return to user
  out <- list(n=n,m=m,d=d,
              codes=codes,
              x=x,
              beta=b,
              xbar=xbar,
              betabar=betabar,
              call=cl)

  class(out) <- c("ideal")

  ## and, finally, if the user wanted meanzero
  if(normalize){
    if(verbose)
      cat("...normalizing output (post-processing)...")
    out <- postProcess(out,
                       constraints="normalize")
    if(verbose)
      cat("done\n")
  }  
  return(out) 
}

x.startvalues <- function(x,d,scale=TRUE,constraint=NULL,verbose=FALSE){
  if(verbose)
    cat("will use eigen-decomposition method to get start values for ideal points...")

  ## from Jong Hee Park
  row.mean <- apply(x, 1, mean, na.rm=TRUE)
  col.mean <- apply(x, 2, mean, na.rm=TRUE)
  dc1 <- sweep(x, 1, row.mean)
  dc2 <- sweep(dc1, 2, col.mean)
  dc <- dc2 + mean(x, na.rm = T)
  
  r <- cor(t(dc),use="pairwise")
  r[is.na(r)] <- 0
  e <- eigen(r)
  v <- e$vectors[,1:d]
  v <- as.matrix(v)
  if(scale){
    for(i in 1:d){
      v[,i] <- v[,i]*sqrt(e$value[i])
    }
  }

  if (!is.null(constraint)) {
    v <- predict(lm(constraint ~ v), newdata=as.data.frame(v)) 
  }
  if(verbose)
    cat("done\n")
  return(v)
}

probit <- function(y,x){
  glmobj <- glm(y ~ x,
                family=binomial(link=probit))
  b <- coef(glmobj)
  k <- length(b)
  b <- b[c(2:k,1)]   ## put intercept last
  b
}

b.startvalues <- function(v,x,d,verbose=FALSE){
  m <- dim(v)[2]
  if(verbose)
    cat(paste("running",
              m,
              "vote-specific probit GLMs\n",
              "for start values for item/bill parameters\n",
              "conditional on start values for ideal points..."))
  
  b <- matrix(NA,m,d+1)
  for(j in 1:m){
    b[j,] <- probit(y=v[,j],x=x)
  }

  ## check for crazy discrimination parameters
  for(j in 1:d){
    bad <- is.na(b[,j])
    b[bad,j] <- 0
  }
  b[,d+1] <- -b[,d+1]     ## flip the sign on the intercepts, make it a difficulty parameter
  if(verbose)
    cat("done\n")
  b
}

