rollcall <- function(data,
                     yea=1,nay=0,
                     missing=NA,
                     notInLegis=9,
                     legis.names=NULL,
                     vote.names=NULL,
                     legis.data=NULL,
                     vote.data=NULL,
                     desc=NULL,
                     source=NULL){
  
  ## codes, check and package
  codes <- list()
  if(!is.null(yea))
    codes$yea <- yea
  else
    codes$yea <- NULL
  if(!is.null(nay))
    codes$nay <- nay
  else
    codes$nay <- NULL
  if(!is.null(notInLegis))
    codes$notInLegis <- notInLegis
  else
    codes$notInLegis <- NULL
  if(!is.null(missing))
    codes$missing <- missing
  else
    codes$missing <- NULL
  if(checkCodes(codes))
    stop("codes are not unique\n")

  ## get a roll call matrix from input
  v <- NULL
  if((is.list(data)) && ("votes" %in% names(data))){ 
    v <- data$votes
  }
  else {
    v <- data
  }
  
  if (!is.matrix(v)) {
    v <- as.matrix(v)
  }
  
  ## check votes
  if(checkVotes(v,codes))
    stop("rollcall: bad votes")
  
  ## identifying tags for legislators
  nm <- legis.names
  if(is.null(nm)) ## look for legis.names var in data
    if(is.list(data))
      if(!is.null(data$legis.names))
        nm <- data$legis.names
  
  if(!is.null(nm)){  ## check that any names found by here are ok
    if(length(unique(nm))!=nrow(v)){
      cat("supplied legislator names do not match number of rows\n")
      cat("in roll call matrix; will use default names\n")
      nm <- NULL
    }
  }
  
  if(is.null(nm)){  ## make names
    nm <- paste("Legislator",1:nrow(v))
  }
  rownames(v) <- nm
  
  ## vote labels
  lbl <- vote.names
  if(is.null(lbl))
    if(is.list(data))
      if(!is.null(data$vote.names))  
        lbl <- data$vote.names

  if(!is.null(lbl)){  ## check that vote names are ok
    if(length(unique(lbl))!=ncol(v)){
      cat("supplied vote names do not match number of columns\n")
      cat("in roll call matrix; will use default names\n")
      lbl <- NULL
    }
  }
  if(is.null(lbl)){  ## make name
    lbl <- paste("Vote",1:ncol(v))
  }
  colnames(v) <- lbl
  
  ## legislator attributes
  if(!is.null(legis.data)){
    if(nrow(legis.data)!=nrow(v))
      stop("legislator data does not match number of legislators in roll call matrix")
    else
      rownames(legis.data) <- nm
  }

  if(!is.null(vote.data)){
    if(nrow(vote.data)!=ncol(v))
      stop("rows in vote.data does not match number of votes in roll call matrix")
    else
      rownames(vote.data) <- lbl
  }

  ## description of roll call voting matrix
  dsc <- desc
  if(is.list(data))
    if(!is.null(data$desc))
      dsc <- data$desc

  ## package up for output
  out <- list(votes=v,
              codes=codes,
              n=dim(v)[1],
              m=dim(v)[2],
              legis.data=legis.data,
              vote.data=vote.data,
              desc=dsc,
              source=source)

  class(out) <- c("rollcall")
  out
}

## loop over roll call matrix
## find which votes to DROP (lop-sided)
lopfunc <- function(x,lop=NULL){
  n <- sum(!is.na(x))
  if(is.null(lop))  ## throw away unanimous votes
    toss <- (sum(x, na.rm=T)==n) || (sum(x,na.rm=T)==0)
  else
    toss <- (sum(x==1,na.rm=T)/n <= lop) || (sum(x==0,na.rm=T)/n <= lop)
  
  toss
}

printDescription <- function(object){
  if(inherits(object,"rollcall") & !is.null(object$desc))
    cat(paste("Description:\t",object$desc,"\n"))
  invisible(NULL)
}

printSource <- function(object){
  if(inherits(object,"rollcall") & !is.null(object$source))
    cat(paste("Source:\t\t",object$source,"\n"))
  invisible(NULL)
}

printHeader <- function(object){
  if(inherits(object,"rollcall")){
    printDescription(object)
    printSource(object)
  }
  invisible(NULL)
}

print.rollcall <- function(x,print.votes=FALSE, ...){
  printHeader(x)
  
  cat(paste("Number of Legislators:\t",x$n,"\n"))
  cat(paste("Number of Votes:\t",x$m,"\n"))

  cat("\n")
  printCodes(x$codes)
  cat("\n")

  if(!is.null(x$legis.data)){
    cat("Legislator-specific variables:\n")
    print(names(x$legis.data))
  }
  if(!is.null(x$vote.data)){
    cat("Vote-specific variables:\n")
    print(names(x$vote.data))
  }
  
  if (print.votes)
   print(x$votes) 

  cat(paste("Detailed information is available via the summary function.\n"))
  
  invisible(NULL)
}

## check Votes
checkVotes <- function(object,codes=object$codes){
  if(inherits(object, "rollcall")){
    mat <- object$votes
  } else {
    if(inherits(object, "matrix")) {
      mat <- object
    }
  }
  if(is.null(codes))
    stop("checkVotes: no codes supplied")
  
  flag <- FALSE
  if(!all(mat[!is.na(mat)] %in% unlist(codes))){
    cat("checkVotes: Your data contains values other than the codes\n")
    cat("checkVotes: you supplied as representing Yea, Nay, Missing or\n")
    cat("checkVotes: not in legislature.\n")
    cat("checkVotes: You specified:\n")
    cat(paste("Yea: ",
              paste(as.character(codes$yea),collapse=" "),
              "\n",
              "Nay: ",
              paste(as.character(codes$nay),collapse=" "),
              "\n",
              "Missing: ",
              paste(as.character(codes$missing),collapse=" "),
              "\n",
              "Not In Legislature: ",
              paste(as.character(codes$notInLegis),collapse=" "),
              "\n",sep=""))
    cat("checkVotes: Your data has the following unique values and frequency counts:\n")
    print(table(mat,exclude=NULL))
    cat("\n")
    flag <- TRUE
  }
  flag
}

## check Codes
checkCodes <- function(codes){
  flag <- FALSE
  n <- length(codes)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      common <- intersect(codes[[i]],codes[[j]])
      if(length(common)!=0)
        flag <- TRUE
    }
  }
  flag
}

## show codes on-screen
printCodes <- function(codes){
  cat("Using the following codes to represent roll call votes:\n")
  if(!is.null(codes$yea))
    cat(paste("Yea:\t\t",
              paste(as.character(codes$yea),collapse=" "),
              "\n"))
  else
    cat(paste("Yea: <no code supplied -- this is almost surely an error>\n"))
  if(!is.null(codes$nay))
    cat(paste("Nay:\t\t",
              paste(as.character(codes$nay),collapse=" "),
              "\n"))
  else
    cat(paste("Nay: <no code supplied -- this is almost surely an error>\n"))
  if(!is.null(codes$missing))
    cat(paste("Abstentions:\t",
              paste(as.character(codes$missing),collapse=" "),
              "\n"))
  else
    cat(paste("Abstentions: <no code supplied -- this is almost surely an error>\n"))
  if(!is.null(codes$notInLegis))
    cat(paste("Not In Legislature:\t",
              paste(as.character(codes$notInLegis),collapse=" "),
              "\n"))
  invisible(NULL)
}

## convert codes to 0 and 1 etc
convertCodes <- function(object,codes=object$codes){
  if(!is.list(codes))
    stop("convertCodes: codes needs to be a list")
  if(is.null(codes$yea) | is.null(codes$nay))
    stop("convertCodes: no Yea and/or Nay code supplied")
  if(!is.matrix(object$votes))
    stop("convertCodes: supplied rollcalls are not in a matrix")

  ## conversions
  theCodes <- codes
  tmp <- matrix(-999,
                dim(object$votes)[1],
                dim(object$votes)[2])
  dimnames(tmp) <- dimnames(object$votes)
  tmp[object$votes %in% theCodes$yea] <- 1
  tmp[object$votes %in% theCodes$nay] <- 0
  if(!is.null(theCodes$missing)){
    if(!any(is.na(theCodes$missing)))
      theCodes$missing <- c(theCodes$missing,NA)
    tmp[object$votes %in% theCodes$missing] <- NA
  }
  if(!is.null(theCodes$notInLegis))
    tmp[object$votes %in% theCodes$notInLegis] <- NA
  
  bad <- tmp[!is.na(tmp)] == -999
  if(any(bad)){
    cat("convertCodes: not all rollcall votes converted to 0, 1, NA.\n")
    cat("convertCodes: information in codes not exhaustive.\n")
    cat(paste("convertCodes: setting remaining",sum(bad),"votes to NA\n"))
    tmp[bad] <- NA
  }
  tmp
}

## for each legislator compute how often they vote with the direction
## in which a majority of their party voted
partyLoyalty <- function(object){
  theParties <- unique(object$legis.data$party)
  nParties <- length(theParties)

  ## what was the majority direction by vote, by party
  majOutcome <- function(x){
    as.numeric(names(which.max(table(x))))
  }
  
  partyDirections <- matrix(NA,object$m,nParties)
  for(p in 1:nParties){
    thisParty <- object$legis.data$party==theParties[p]
    ## only do this if the party is bigger than two legislators
    if(sum(thisParty,na.rm=TRUE)>2){
      foo <- apply(object$votes[thisParty,],2,
                   majOutcome)
      if(is.list(foo)){
        foo[which(lapply(foo,length)==0)] <- NA
        foo <- unlist(foo)
      }
      partyDirections[,p] <- foo
     }
   }
  dimnames(partyDirections) <- list(dimnames(object$votes)[[2]],
                                    theParties)
  
  ## now compare individual voting histories
  ## with voting scores
  partyLoyalty <- rep(NA,object$n)
  legisParty <- match(x=object$legis.data$party,
                      table=theParties)
  
  goodCompare <- function(x,y){
    ok <- !is.na(x) & !is.na(y)
    goodmatch <- sum(x[ok] == y[ok]) + sum(is.na(x) & is.na(y))
    out <- goodmatch/length(x) * 100
    out
  }
  
  for(i in 1:object$n){   ## loop over legislators
    ## dont do it where partyDirections are undefined
    ## (i.e., small numbers of indeps etc)
    if(!all(is.na(partyDirections[,legisParty[i]]))){  
      theDirections <- partyDirections[,legisParty[i]]  
      partyLoyalty[i] <- goodCompare(theDirections,
                                     object$votes[i,])
    }
  }
  partyLoyalty
}

lopLook <- function(margins,cutOff){
  extremeMat <- rep(NA,cutOff+1)
  for(j in 0:cutOff){
    extremeMat[j+1] <- sum(margins[,1]==j | margins[,2]==j)
  }
  extremeMat
}

vectorRepresentation <- function(object,
                                 dropList=list(codes=c("missing",
                                                 "notInLegis"))){
  if(!inherits(object, "rollcall"))
    stop("vectorRepresentation only defined for objects of class rollcall")
  if(is.null(object$codes))
    stop("no rollcall codes")
  else
    codes <- object$codes
  
  if(is.null(dropList) |
     length(dropList)==0 |
     is.null(dropList$codes) |
     length(dropList$codes)==0){
    cat("missing arguments for drop, vectorRepresentation will use defaults")
    dL <- list(codes=c("missing","notInLegis"))
  }
  else
    dL <- dropList
  tmpRollCall <- dropRollCall(object,dropList=dL)

  badCodes <- match(dL$codes,names(codes))
  if(any(is.na(badCodes)))
    stop("couldn't find codes to drop\n")
  else{
    dropCodes <- NULL
    for(j in badCodes)
      dropCodes <- c(dropCodes,codes[[j]])
  }
  cat(paste("vectorRepresentation: dropCodes=",
            paste(dropCodes,collapse=", "),
            "\n"))
  
  n <- tmpRollCall$n
  m <- tmpRollCall$m
  v <- tmpRollCall$votes
  y <- matrix(NA,n*m,3)
  dimnames(y)[[2]] <- list("vote","i","j")

  z <- 1
  for(i in 1:n){
    for(j in 1:m){
      vij <- v[i,j]
      if(!(vij %in% dropCodes)){
        if(vij %in% codes$yea)
          y[z,1] <- 1
        if(vij %in% codes$nay)
          y[z,1] <- 0
        y[z,2] <- i
        y[z,3] <- j
        z <- z + 1
      }
    }
  }
  y <- y[1:(z-1),]
  y
}

matchDimnames <- function(labs,codes){
  codesNames <- names(codes)
  nLabs <- length(labs)
  whichCode <- rep(NA,nLabs)
  for(i in 1:nLabs){
    tmp <- unlist(lapply(codes,function(x)as.numeric(labs[i]) %in% x))
    #print(tmp)
    if(sum(tmp)==1)
      whichCode[i] <- which.max(tmp)
  }
  out <- codesNames[whichCode]
  out <- paste(labs," (",out,")",sep="")
  out
}

summary.rollcall <- function(object,
                             dropList=NULL,
                             ##list(codes="notInLegis",
                             ##  lop=0),
                             verbose=FALSE,
                             debug=FALSE,
                             ...){
  if(!inherits(object, "rollcall"))
    stop("summary.rollcall only operates on objects of class rollcall")

  mc <- match.call()   ## how were we called
  if(is.null(mc$dropList))
    mc$dropList <- dropList
  if(is.null(mc$verbose))
    mc$verbose <- verbose
  
  legisTab <- NULL
  voteTab <- NULL
  dropTab <- NULL
  partyLoyaltyScores <- NULL
  lopSided <- NULL

  if(!is.null(object$dropInfo)){
    cat(paste("The input rollcall object already has a dropInfo component,\n",
              "meaning that it is the product of dropRollCall.\n",
              "This summary and the execution of any current dropList\n",
              "proceeds conditional on the previous dropList.\n"),
        sep="")
  }   

  ## process user options re dropping votes/legislators etc
  if(!is.null(dropList)){
    tmpRollCall <- dropRollCall(object,dropList,debug=debug)
  }
  else
    tmpRollCall <- object

  v <- tmpRollCall$votes
  
  ## party breakdown, if available
  haveParty <- !is.null(tmpRollCall$legis.data$party)
  partyTab <- NULL
  if(haveParty)
    partyTab <- table(tmpRollCall$legis.data$party,exclude=NULL)

  ## get any exclude codes
  
  allVotes <- table(v)
  allVotes <- cbind(allVotes,
                    allVotes/sum(allVotes)*100)
  dimnames(allVotes)[[2]] <- c("Count","Percent")
  dimnames(allVotes)[[1]] <- matchDimnames(dimnames(allVotes)[[1]],
                                           tmpRollCall$codes)
  
  ## what was clobbered by dropRollCall
  if(!is.null(tmpRollCall$dropInfo))
    if("new" %in% names(tmpRollCall$dropInfo))
      dropTab <- lapply(tmpRollCall$dropInfo$new[c("legislators","votes")],
                        table,exclude=NULL)
    else
      dropTab <- lapply(tmpRollCall$dropInfo[c("legislators","votes")],
                        table,exclude=NULL)
  
  if(verbose){
    ## breakdowns by legislator
    if(debug)
      cat("computing breakdowns by legislator...")
    legisTab <- t(apply(v,1,
                        marginWithCodes,
                        codes=tmpRollCall$codes))
    dimnames(legisTab)[[2]] <- c(names(tmpRollCall$codes),
                                 "Total",
                                 paste(names(tmpRollCall$codes),"%",sep=""))
    ## breakdowns by vote
    cat("by vote...")
    voteTab <- t(apply(v,2,
                       marginWithCodes,
                       codes=tmpRollCall$codes))
    dimnames(voteTab)[[2]] <- dimnames(legisTab)[[2]]
    
    lopSided <- lopLook(voteTab,floor(.05*tmpRollCall$n))
    names(lopSided) <- as.character(0:floor(.05*tmpRollCall$n))

    ## party loyalty
    if(haveParty){
      cat("and party loyalty scores")
      partyLoyaltyScores <- partyLoyalty(tmpRollCall)
    }
    cat("\n")
  }
  
  out <- list(n=tmpRollCall$n,
              m=tmpRollCall$m,
              codes=tmpRollCall$codes,
              allVotes=allVotes,
              partyTab=partyTab,
              lopSided=lopSided,
              legisTab=legisTab,
              dropTab=dropTab,
              partyLoyalty=partyLoyaltyScores,
              voteTab=voteTab,
              call=mc)
  class(out) <- "summary.rollcall"
  out
}

printDropTab <- function(x){
  for(i in 1:length(x))
    if("FALSE" %in% names(x[[i]])){
      cat(paste("dropRollCall deleted",
                x[[i]]["FALSE"],
                "of",
                sum(x[[i]]),
                names(x)[[i]],
                "\n"))
    }
    else{
      cat(paste("dropRollCall deleted no",names(x)[[i]],"\n"))
    }
  invisible(NULL)
}

print.summary.rollcall <- function(x, digits=1, ...){
  if(!inherits(x, "summary.rollcall"))
    stop("print.summary.rollcall only defined for objects of class summary.rollcall")
    
  rcObj <- x$call$object
  verbose <- x$call$verbose
  if(is.null(eval(rcObj)))
    stop("can't find rollcall object")

  if(length(rcObj)>1)
    rcObjName <- format(rcObj)
  else
    rcObjName <- rcObj
  if(!is.null(eval(rcObj)$desc))
    cat(paste("\nSummary of rollcall object",
              rcObjName,
              "\n\n"))
  printHeader(eval(rcObj))

  cat(paste("\nNumber of Legislators:\t\t",x$n))
  cat(paste("\nNumber of Roll Call Votes:\t",x$m))
  cat("\n\n")

  ## was summary called from dropRollCall directly? 
  ## if so then dump the drop info from the object 
  
  if(!is.null(x$dropTab))
    printDropTab(x$dropTab)

  ## if(!is.null(x$call$dropList))
##     cat(paste("This summary ignores voting decisions that are coded ",
##               paste(x$call$dropList$codes,collapse=" or "),
##               ".\n",sep=""))
##   if(!is.null(x$call$dropList$lop))
##     if(x$call$dropList$lop==0)
##       cat("This summary computed after dropping unanimous roll calls.\n")
##     else
##       cat(paste("This summary computed after dropping roll calls with",
##                 x$call$dropList$lop,
##                 "or fewer legislators voting\nin the minority.\n"))
##   if(length(x$call$dropList)>2){
##     cat("Other restrictions are being applied. The full dropList is:\n")
##     print(x$call$dropList)
##   }
     
  cat("\n")
  printCodes(x$codes)
  cat("\n")
  
  if(!is.null(x$partyTab)){
    cat("Party Composition:")
    print(x$partyTab)
  }
  cat("\nVote Summary:\n")
  print(round(x$allVotes,1))

  if(!is.null(x$lopSided) | !all(x$lopSided==0)){
    cat("\nLop-sided Votes (Number Voting in Minority), and Frequencies:\n")
    print(x$lopSided)
  }

  if(!verbose)
    cat(paste("\nUse summary(",
              rcObjName,
              ",verbose=TRUE) for more detailed information.\n",sep=""))
  
  if(verbose){
    if(!is.null(x$partyLoyalty)){
      cat("\nSummary By Legislator: Counts, Percentages and Party Loyalty\n")
      foo <- cbind(round(x$legisTab,digits),
                   round(x$partyLoyalty))
      dimnames(foo)[[2]][ncol(foo)] <- "Party Loyalty"
    }
    else{
      cat("\nSummary By Legislator: Counts and Percentages\n")
      foo <- round(x$legisTab,digits)
    }
    print(foo)
    
    cat("\nSummary By Vote\n")
    print(round(x$voteTab,digits))
    cat("\n")
  }

  invisible(NULL)
}

