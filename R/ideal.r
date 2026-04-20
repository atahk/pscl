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
                  use.voter=NULL,
                  n.chains=1L,
                  mc.cores=getOption("mc.cores", 1L),
                  seed=NULL,
                  chain.sd=1.0,
                  n.threads=getOption("pscl.n.threads", 1L),
                  .backend=c("classic","v2","v3","v4","v5")){

  .backend <- match.arg(.backend)
  .c_entry <- switch(.backend,
                     classic = "IDEAL",
                     v2      = "IDEAL_v2",
                     v3      = "IDEAL_v3",
                     v4      = "IDEAL_v4",
                     v5      = "IDEAL_v5")
  n.threads <- max(1L, as.integer(n.threads))
  if (.backend != "v5" && n.threads > 1L) {
    warning("n.threads > 1 is only used by .backend='v5'; ignored.\n")
    n.threads <- 1L
  }

  n.chains <- as.integer(n.chains)
  if (is.na(n.chains) || n.chains < 1L)
    stop("n.chains must be a positive integer\n")
  mc.cores <- max(1L, as.integer(mc.cores))

  ## Oversubscription check: warn if (concurrent chains) * (threads per chain)
  ## exceeds the number of physical cores. We don't auto-budget -- the user
  ## may be running on a shared machine, or alongside other jobs. Stan and
  ## cmdstanr take the same hands-off stance.
  ncores <- tryCatch(parallel::detectCores(logical = FALSE),
                     error = function(e) NA_integer_)
  if (!is.na(ncores)) {
    workers_eff <- min(mc.cores, n.chains)
    demand <- workers_eff * n.threads
    if (demand > ncores) {
      warning(sprintf(
        "min(mc.cores, n.chains) * n.threads = %d exceeds %d physical core(s); consider lowering one of them to avoid oversubscription.\n",
        demand, ncores))
    }
  }

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
  if(is.null(cl$n.chains))
    cl$n.chains <- n.chains

  mda <- FALSE

  ## check validity of user arguments
  if (!inherits(object, "rollcall")){
    stop("object must be of class rollcall\n")
  }
  if(((d%%1) != 0) || (d<1)){
    stop("d is not a positive integer\n")
  }

  if(((thin%%1)!=0) || (thin<1)) {
    stop("thin is not a positive integer\n")
  }

  if(((maxiter%%1)!=0) || (maxiter<1)) {
    stop("maxiter is not a positive integer\n")
  }

  if(!is.list(dropList)){
    stop("dropList must be a list\n")
  }

  if(!is.list(codes)){
    stop("codes must be a list\n")
  }

  ##check iterations and thinning
  if ((maxiter%%thin)!=0) {
    stop("maxiter must be a multiple of thin\n")
  }

  if ((burnin%%thin)!=0) {
    stop("burnin must be a multiple of thin\n")
  }

  if (burnin >= maxiter){
    stop("burnin must be less than maxiter\n")
  }

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
    if(checkCodes(tmpObject$codes)){
      stop("supplied codes fail redundancy checks\n")
    }
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
      if (!is.vector(use.voter)){
          stop("use.voter must be a vector of length n\n")
      }
      if (n != length(use.voter)){
          stop("use.voter must be a vector of length n\n")
      }
  }

  ## map roll call votes into binary format required by ideal
  if(verbose){
    printCodes(codes)
    cat("\n")
  }

  if(checkVotes(y$votes,codes)){
    stop("rollcall: can't map all votes to 0|1 using supplied codes\n")
  }

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

  ## ############################################
  ## per-chain seeds and perturbed start values
  ## ############################################
  if (n.chains > 1L) {
    if (is.null(seed)) {
      seeds <- sample.int(.Machine$integer.max, n.chains)
    } else if (length(seed) == 1L) {
      set.seed(as.integer(seed))
      seeds <- sample.int(.Machine$integer.max, n.chains)
    } else {
      if (length(seed) != n.chains)
        stop("seed must be NULL, a single integer, or a vector of length n.chains\n")
      seeds <- as.integer(seed)
    }

    chain_starts <- lapply(seq_len(n.chains), function(i) {
      if (i == 1L) {
        list(x = xstart, b = bstart)
      } else {
        list(x = xstart + rnorm(length(xstart), sd = chain.sd),
             b = bstart + rnorm(length(bstart), sd = chain.sd))
      }
    })
  } else {
    if (!is.null(seed)) set.seed(as.integer(seed))
    seeds <- NA_integer_
    chain_starts <- list(list(x = xstart, b = bstart))
  }

  ## ############################################
  ## helper: run one chain through .C and assemble the ideal object
  ## ############################################
  runOneChain <- function(chain_idx) {
    if (n.chains > 1L) set.seed(seeds[chain_idx])
    cs <- chain_starts[[chain_idx]]
    xs <- cs$x
    bs <- cs$b

    file_i <- file
    if (usefile && n.chains > 1L) {
      file_i <- sub("(\\.[^./\\\\]+)?$",
                    sprintf("_chain%d\\1", chain_idx),
                    file)
    }

    chain_verbose <- verbose && n.chains == 1L

    if (n.chains == 1L)
      cat("\nStarting MCMC Iterations...\n")

    if (usefile) {
      if (length(legis.names) == n) {
        cat(paste("\"", c("Iteration", legis.names), "\"",
                  sep="", collapse=","),
            file=file_i)
      } else {
        cat(paste("\"", c("Iteration", paste("x", 1:n, sep="")), "\"",
                  sep="", collapse=","),
            file=file_i)
      }
      if (store.item){
        cat(",",
            paste("\"",
                  c(paste("b",
                          as.vector(apply(expand.grid(1:m,1:(d+1)),1,paste,collapse=".")),
                          sep=".")),"\"",
                  sep="", collapse=","),
            sep="", file=file_i, append=TRUE)
      }
      cat("\n", file=file_i, append=TRUE)

      .C_args <- list(.c_entry,
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
                      as.double(xs),
                      as.double(bs),
                      xoutput=0,
                      boutput=0,
                      as.integer(burnin),
                      usefile,
                      as.logical(store.item),
                      as.character(file_i),
                      as.logical(chain_verbose),
                      as.logical(!is.null(use.voter)),
                      as.integer(use.voter))
      if (.backend == "v5") .C_args <- c(.C_args, list(as.integer(n.threads)))
      output <- do.call(.C, .C_args)
    }
    else if (!store.item) {
      .C_args <- list(.c_entry,
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
                      as.double(xs),
                      as.double(bs),
                      xoutput=as.double(rep(0,n*d*numrec)),
                      boutput=0,
                      as.integer(burnin),
                      usefile,
                      as.logical(store.item),
                      as.character(file_i),
                      as.logical(chain_verbose),
                      as.logical(!is.null(use.voter)),
                      as.integer(use.voter))
      if (.backend == "v5") .C_args <- c(.C_args, list(as.integer(n.threads)))
      output <- do.call(.C, .C_args)
    }
    else {
      .C_args <- list(.c_entry,
                      PACKAGE=.package.Name,
                      as.integer(n), as.integer(m), as.integer(d), as.double(yToC),
                      as.integer(maxiter), as.integer(thin), as.integer(impute),
                      as.integer(mda),
                      as.double(xp), as.double(xpv), as.double(bp),
                      as.double(bpv), as.double(xs), as.double(bs),
                      xoutput=as.double(rep(0,n*d*numrec)),
                      boutput=as.double(rep(0,m*(d+1)*numrec)),
                      as.integer(burnin),
                      usefile,
                      as.logical(store.item),
                      as.character(file_i),
                      as.logical(chain_verbose),
                      as.logical(!is.null(use.voter)),
                      as.integer(use.voter))
      if (.backend == "v5") .C_args <- c(.C_args, list(as.integer(n.threads)))
      output <- do.call(.C, .C_args)
    }

    if (n.chains == 1L) cat("\n")

    ## parse returns from C
    xbar <- NULL
    betabar <- NULL
    if (!usefile) {
      itervec <- seq(burnin,maxiter,by=thin)
      keep <- itervec > burnin

      xx <- array(output$xoutput, c(n,d,numrec))
      xx <- aperm(xx, c(3,1,2))
      dimnames(xx) <- list(itervec,
                           legis.names,
                           paste("D",1:d,sep=""))
      if(chain_verbose)
        cat("...computing posterior means for ideal points...")
      xbar <- getMean(keep, xx)
      if(chain_verbose)
        cat("done\n")

      if(store.item){
        bb <- array(output$boutput, c(m,d+1,numrec))
        dimnames(bb) <- list(vote.names,
                             c(paste("Discrimination D",1:d,sep=""),
                               "Difficulty"),
                             itervec)
        bb <- aperm(bb, c(3,1,2))
        if(chain_verbose)
          cat("...computing posterior means for item parameters...")
        betabar <- getMean(keep, bb)
        if(chain_verbose)
          cat("done\n")
      } else {
        bb <- NULL
      }
    } else {
      bb <- xx <- NULL
    }

    out <- list(n=n, m=m, d=d,
                codes=codes,
                x=xx,
                beta=bb,
                xbar=xbar,
                betabar=betabar,
                call=cl)
    class(out) <- "ideal"

    if (normalize) {
      if (chain_verbose)
        cat("...normalizing output (post-processing)...")
      out <- postProcess(out, constraints="normalize")
      if (chain_verbose)
        cat("done\n")
    }
    out
  }

  ## ############################################
  ## dispatch: single chain vs multi-chain
  ## ############################################
  if (n.chains == 1L) {
    return(runOneChain(1L))
  }

  ## helper for the multi-chain auto-align fallback below
  flipOne <- function(f) {
    f$x    <- -f$x
    f$xbar <- -f$xbar
    if (!is.null(f$beta)) {
      f$beta[,,1]    <- -f$beta[,,1]
      f$betabar[,1]  <- -f$betabar[,1]
    }
    f
  }

  workers <- min(mc.cores, n.chains)
  use_parallel <- workers > 1L
  use_fork <- use_parallel && .Platform$OS.type != "windows"

  cat(sprintf("Running %d chains", n.chains))
  if (use_parallel)
    cat(sprintf(" in parallel on %d %s worker%s",
                workers,
                if (use_fork) "fork" else "PSOCK",
                if (workers > 1L) "s" else ""))
  cat("...\n")

  if (use_fork) {
    fits <- parallel::mclapply(seq_len(n.chains),
                               runOneChain,
                               mc.cores = workers,
                               mc.preschedule = FALSE)
    bad <- vapply(fits, inherits, logical(1), what = "try-error")
    if (any(bad))
      stop("chain(s) failed: ", paste(which(bad), collapse=", "))
  } else if (use_parallel) {
    ## PSOCK cluster (Windows, or user-requested portable path).
    ## parLapply() serialises runOneChain together with its enclosing
    ## environment, so all free variables (priors, starts, yToC, seeds,
    ## ...) travel to each worker. The worker needs pscl loaded so that
    ## .package.Name and the native entry points resolve.
    cl_ <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl_), add = TRUE)
    parallel::clusterEvalQ(cl_, {
      suppressPackageStartupMessages(library(pscl))
      NULL
    })
    fits <- parallel::parLapply(cl_, seq_len(n.chains), runOneChain)
  } else {
    fits <- lapply(seq_len(n.chains), runOneChain)
  }
  cat(sprintf("All %d chains complete.\n", n.chains))

  ## chain-to-chain auto-alignment at d=1: if the user did not identify
  ## polarity via constrain.legis()-supplied priors, chains can land in
  ## mirror-image modes. Use cor(xbar) vs chain 1 as a cheap fallback so
  ## downstream idealToMCMC() / gelman.diag() aren't artificially inflated.
  ## With identifying priors, chains are already aligned and cor > 0 so
  ## nothing flips -- safe to always run.
  if (d == 1L && length(fits) > 1L) {
    x1 <- as.vector(fits[[1]]$xbar)
    for (i in 2:length(fits)) {
      if (cor(x1, as.vector(fits[[i]]$xbar)) < 0)
        fits[[i]] <- flipOne(fits[[i]])
    }
  }

  class(fits) <- c("idealList", "list")
  attr(fits, "seeds")    <- seeds
  attr(fits, "n.chains") <- n.chains
  attr(fits, "call")     <- cl
  fits
}

## ##################################################################
## methods for class idealList
## ##################################################################
print.idealList <- function(x, ...) {
  cat(sprintf("List of %d ideal() chains\n", length(x)))
  if (!is.null(x[[1]]$n))
    cat(sprintf("  %d legislators x %d items x d=%d\n",
                x[[1]]$n, x[[1]]$m, x[[1]]$d))
  seeds <- attr(x, "seeds")
  if (!is.null(seeds))
    cat(sprintf("  seeds: %s\n", paste(seeds, collapse=", ")))
  cat("\nUse idealToMCMC() to obtain a coda mcmc.list for convergence diagnostics.\n")
  invisible(x)
}

summary.idealList <- function(object, ...) {
  lapply(seq_along(object), function(i) {
    cat(sprintf("\n---- chain %d ----\n", i))
    summary(object[[i]], ...)
  })
  invisible(object)
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
