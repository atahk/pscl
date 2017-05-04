## parse and execute drop list used by summary.rollcall and ideal
printDropList <- function(list){
  cat("Dropping elements of rollcall matrix using the following dropList:\n")
  if(!is.null(list$codes)){
    cat("  Voting decisions with the following codes will be set to NA:\n")
    print(list$codes)
    cat("\n")
  }

  if(!is.null(list$lop) & is.numeric(list$lop)){
    if(list$lop==0)
      cat("  Unanimous votes will be dropped.\n")
  else
    cat("  Votes with",list$lop,
        "or fewer legislators voting in the minority",
        "will be dropped.\n")
  }

  if(!is.null(list$legisMin) & is.numeric(list$legisMin))
    cat("  Legislators with",list$legisMin,
        "or fewer non-missing voting decisions",
        "will be dropped.\n")

  if(!is.null(list$dropLegis))
    cat(paste("  Legislators for whom the condition\n",
              "  ",
              deparse(list$dropLegis),
              "\n",
              "  is true (evaluated in the legis.data data frame) will be dropped.\n",
              sep=""))

  if(!is.null(list$dropVotes))
    cat(paste("  Votes for which the condition\n",
              deparse(list$dropVotes),
              "\n",
              "  is true (evaluated in the vote.data data frame) will be dropped.\n",
              sep=""))
  cat("\n")

  invisible(NULL)
}

compareRollCallObjects <- function(old,new){
  legislators <- dimnames(old$votes)[[1]]
  n <- length(legislators)
  votes <- dimnames(old$votes)[[2]]
  m <- length(votes)
  
  if(all(dim(old$votes)==dim(new$votes)))
    keep <- list(legislators=rep(TRUE,n),votes=rep(TRUE,m))
  else{
    newLegis <- dimnames(new$votes)[[1]]
    newVotes <- dimnames(new$votes)[[2]]
    keep <- list(legislators=legislators%in%newLegis,
                 votes=votes%in%newVotes)
  }
  names(keep$legislators) <- legislators
  names(keep$votes) <- votes

  keep 
}


dropRollCall <- function(object,dropList=NULL,debug=FALSE){
  if(class(object)!="rollcall"){
    stop("dropRollCall only works for objects of class rollcall.")
  }

  tmpRollCall <- object
  
  if(!is.list(dropList) | is.null(dropList)){
    cat("dropList must be a non-null list or alist.\nNo subsetting will occur.\n")
    return(object)
  }

  if(debug)
    printDropList(dropList)

  flag <- TRUE   ## raise the flag
  counter <- 1
  while(flag){   ## loop until the flag goes down
    if(debug)
      cat(paste("\ndropRollCall: Pass number",counter,"over roll call object\n"))

    v <- tmpRollCall$votes
    dimOld <- dim(v)               ## store this
    if(debug)
      cat(paste("  The roll call matrix has dimension",
                dimOld[1],"legislators",
                "and",
                dimOld[2],"rollcalls.\n"))

    ## strip out user-designated votes of a particular code
    if(!is.null(dropList$codes) & length(dropList$codes)>0){
      if(debug)
        cat("  Processing dropList voting codes...\n")
      dc <- dropList$codes
      dCodes <- NULL
      if(all(is.character(dc))){   ## named element of codes list?
        dropCodes <- match(dc,names(tmpRollCall$codes))
        dropCodes <- dropCodes[!is.na(dropCodes)]
        
        if(length(dropCodes)>0){
          for(j in dropCodes)
            dCodes <- c(dCodes,tmpRollCall$codes[j])  ## drop these
          keepCodes <- !(names(tmpRollCall$codes) %in% dc)
          keepCodes <- tmpRollCall$codes[keepCodes]
          tmpRollCall$codes <- keepCodes    
        }
      }
      if(is.numeric(dc)){    ## or numeric elements
        dCodes <- dc[dc %in% unique(as.vector(v))]
      }
      
      bad <- v %in% dCodes
      if(debug)
        cat(paste("  dropRollCall will set",sum(bad),"voting decisions to NA.\n"))
      tmpRollCall$votes[bad] <- NA
      rm(bad)
    }

    dropLegis <- rep(FALSE,dim(v)[1])
    dropVotes <- rep(FALSE,dim(v)[2])
    
    ## drop legislators if too little data
    if(!is.null(dropList$legisMin)){
      if(debug)
        cat("  dropRollCall processing minimum votes by legislator (legisMin) restrictions...\n")
      legisMin <- dropList$legisMin
      if(length(legisMin)!=1 |
         is.na(legisMin) |
         !is.numeric(legisMin) |
         legisMin >= tmpRollCall$m)
        stop("  Bad value for legisMin in drop list.")
      vtmp <- convertCodes(tmpRollCall)
      goodCount <- apply(vtmp,1,function(x)sum(!is.na(x)))
      dropLegis <- dropLegis | goodCount<legisMin
    }
    
    ## check for subsetting in legis.data
    if(!is.null(dropList$dropLegis)){
      r <- dropRollCallViaData(dropList$dropLegis,
                               object=tmpRollCall,
                               d=expression(legis.data))
      if(!is.null(r))
        dropLegis <- dropLegis | r
    }
    
    ## drop votes by lop-sidedness
    if(!is.null(dropList$lop) & is.numeric(dropList$lop)){
      if(debug)
        cat("  dropRollCall processing lop-sided restrictions...\n")
      lop <- dropList$lop
      if(length(lop)!=1 |
         is.na(lop) |
         !is.numeric(lop) |
         lop < 0 |
         lop >= tmpRollCall$n)
        stop("  Invalid value for lop")
      if(is.null(tmpRollCall$voteMargins)){
        if(debug)
          cat("    Computing vote margins...\n")
        tmpRollCall <- computeMargins(tmpRollCall,dropList=NULL)
      }
      r <- tmpRollCall$voteMargins[,"Min"] <= lop
      if(debug)
        cat(paste("    dropRollCall will drop",sum(r),"roll calls due to lop-sidedness.\n"))
      dropVotes <- dropVotes | r
      if(debug)
        cat("  dropRollCall finished processing lop-sided restrictions.\n")
    }
    
    ## check for subsetting in vote.data
    if(!is.null(dropList$dropVotes) & counter==1){
      r <- dropRollCallViaData(dropList$dropVotes,
                               object=tmpRollCall,
                               d=expression(vote.data))
      if(!is.null(r))
        dropVotes <- dropVotes | r
    }

    ## final processing
    if(debug)
      cat(paste("dropRollCall will drop ",
                sum(dropLegis),
                " legislators & ",
                sum(dropVotes),
                " rollcalls.\n",
                sep=""))
    ## if(sum(dropLegis)>0){
##       cat("Dropped Legislators:\n")
##       print(dimnames(tmpRollCall$votes)[[1]][dropLegis])
##     }
##     if(sum(dropVotes)>0){
##       cat("Dropped Votes:\n")
##       print(dimnames(tmpRollCall$votes)[[2]][dropVotes])
##     }

    if(sum(dropLegis>0) | sum(dropVotes)>0)
      tmpRollCall$votes <- tmpRollCall$votes[!dropLegis,!dropVotes]
    
    if(!is.null(tmpRollCall$legis.data) & sum(dropLegis)>0){
      tmpRollCall$legis.data <- tmpRollCall$legis.data[!dropLegis,]
      if(class(object$legis.data)=="data.frame"){
        class(tmpRollCall$legis.data) <- "data.frame"
        names(tmpRollCall$legis.data) <- names(object$legis.data)
      }
    }

    if(!is.null(tmpRollCall$vote.data) & sum(dropVotes)>0){
      tmpRollCall$vote.data <- tmpRollCall$vote.data[!dropVotes,]
      if(class(object$vote.data)=="data.frame"){
        class(tmpRollCall$vote.data) <- "data.frame"
        names(tmpRollCall$vote.data) <- names(object$vote.data)
      }
    }

    if(!is.null(tmpRollCall$voteMargins) & sum(dropVotes)>0)
      tmpRollCall$voteMargins <- tmpRollCall$voteMargins[!dropVotes,]

    dimNew <- dim(tmpRollCall$votes)
    tmpRollCall$n <- dimNew[1]
    tmpRollCall$m <- dimNew[2]

    if(all(dimNew==dimOld)){          ## if no change from previous pass
      if(debug)
        cat("\ndropRollCall has finished processing the rollcall object.\n\n")
      flag <- FALSE                  ## lower the flag, quit loop
    }

    counter <- counter + 1
  }

  ## and finally, add dropped information to rollcall object
  if(!is.null(dropList)){
    newdropInfo <- compareRollCallObjects(object,tmpRollCall)

    if(is.null(tmpRollCall$dropInfo)){
      tmpRollCall$dropInfo <- newdropInfo
      tmpRollCall$dropInfo$dropList <- dropList
    }
    else{
      ## add this dropList to the others
      tmpRollCall$dropInfo <- list(previous=tmpRollCall$dropInfo,
                                   new=list(legislators=newdropInfo$legislators,
                                     votes=newdropInfo$votes,
                                     dropList=dropList))
    }
  }

  return(tmpRollCall)               ## return rollcall object
}

dropRollCallViaData <- function(expr,object,d){
  cf <- match.call()

  f <- try(eval(d,envir=object),silent=TRUE)
  if(inherits(f,"try-error")){
    cat(paste("The data frame ",
              cf$d,
              " was not found in ",
              cf$object,
              ".\n",sep=""))
    cat("Proceeding by ignoring this subsetting restriction.\n")
    return(NULL)
  }

  r <- try(eval(expr,f),silent=TRUE)
  if(inherits(r,"try-error")){
    r <- rep(FALSE,dim(f)[1])
    cat(paste("The assertion ",
              deparse(expr),
              " could not be evaluated in the ",
              cf$d,
              " component of ",
              cf$object,
              ".\n",
              sep=""))
    cat("Proceeding by ignoring this assertion.\n")
  }
  if(!is.logical(r))
    stop("'x' must evaluate to logical")
  r <- r & !is.na(r)
  r
}
