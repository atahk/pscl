## read Keith Poole and Howard Rosenthal formatted data
readKH <- function(file,
                   dtl=NULL,
                   yea=c(1,2,3),
                   nay=c(4,5,6),
                   missing=c(7,8,9),
                   notInLegis=0,
                   desc=NULL,
                   debug=FALSE){

  cat("Attempting to read file in Keith Poole/Howard Rosenthal (KH) format.\n")
  warnLevel <- options()$warn
  options(warn=-1)
  data <- try(readLines(con=file),silent=TRUE)
  if(inherits(data,"try-error")){
    cat(paste("Could not read",file,"\n"))
    return(invisible(NULL))
  }
  options(warn=warnLevel)
  
  cat("Attempting to create roll call object\n")
  voteData <- substring(data,37)
  n <- length(voteData)
  m <- nchar(voteData)[1]
  rollCallMatrix <- matrix(NA,n,m)
  for(i in 1:n){
    rollCallMatrix[i,] <- as.numeric(unlist(strsplit(voteData[i],
                                                     split=character(0))))
  }
  rm(voteData)

  if(!is.null(desc))
    cat(paste(desc,"\n"))
  cat(paste(n,"legislators and",m,"roll calls\n"))
  cat("Frequency counts for vote types:\n")
  tab <- table(rollCallMatrix,exclude=NULL)
  print(tab)

  ## unique numeric identifier for each legislator
  icpsrLegis <- as.numeric(substring(data,4,8))
  
  ## party affiliation
  party <- as.numeric(substring(data,21,23))
  ## convert party to label
  partyfunc <- function(x){
    ##data(partycodes)
    party <- partycodes$party[match(x,partycodes$code)]
    party[party=="Democrat"] <- "D"
    party[party=="Republican"] <- "R"
    party[party=="Independent"] <- "Indep"
    party
  }
  partyName <- partyfunc(party)

  ## convert state ICPSR code to abbreviation
  statename <- function(x){
    ##data(state.info)
    state.info$state[match(x,state.info$icpsr)]
  }
  state <- as.numeric(substring(data,9,10))  ## icpsr code
  KHstateName <- substring(data,13,20)
  stateName <- statename(state)  ## covert to name
  stateAbb <- datasets::state.abb[match(stateName,datasets::state.name)]  ## convert to abbrev
  stateAbb[grep(KHstateName,pattern="^USA")] <- "USA"  ## for presidents

  cd <- as.numeric(substring(data,11,12))
  cdChar <- as.character(cd)
  cdChar[cd==0] <- ""

  ## process legislator names
  lnames <- substring(data,26,36)
  for(i in 1:n){
    lnames[i] <- strip.trailing.space(lnames[i])
    lnames[i] <- strip.after.comma(lnames[i])
  }

  ## finally, produce a tag for each legislator
  legisId <- paste(lnames," (",partyName," ",stateAbb,"-",cdChar,")",sep="")
  legisId <- gsub(x=legisId,pattern="-)",replacement=")")

  ## final check for dups
  ## if we find any, pad with icpsrLegis tag
  if(any(duplicated(legisId))){
    dups <- duplicated(legisId)
    legisId[dups] <- paste(legisId[dups],
                           icpsrLegis[dups])
  }

  ## write legis data
  legis.data <- data.frame(state=stateAbb,
                           icpsrState=state,
                           cd=cd,
                           icpsrLegis=icpsrLegis,
                           party=partyName,
                           partyCode=party)
  dimnames(legis.data)[[1]] <- legisId

  ## do we have a dtl file to read?
  vote.data <- NULL
  if(!is.null(dtl)){
    vote.data <- dtlParser(dtl,debug=debug)
  }
  
  ## finally, call rollcall to assemble working object
  rc <- rollcall(data=rollCallMatrix,
                 yea=yea,
                 nay=nay,
                 missing=missing,
                 notInLegis=notInLegis,
                 legis.names=legisId,
                 legis.data=legis.data,
                 vote.data=vote.data,
                 desc=desc,
                 source=file)
  rc
}

## utility functions
strip.after.comma <- function(x){
  indx <- regexpr(",",x)
  if (indx > 0)
    z <- substring(x,1,indx-1)
  else
    z <- x
  z
}

strip.trailing.space <- function(x){
  indx <- regexpr(" ",x)
  if (indx > 0)
    z <- substring(x,1,indx-1)
  else
    z <- x
  z
}

## read from file, possible web
## readFromFunc <- function(file,debug=TRUE){
##   ## check if this is a URL, starting with such as http, https, or ftp
##   urlStrings <- c("http","ftp")
##   netFile <- any(!is.na(pmatch(urlStrings,file)))
##   if(netFile){
##     if(debug)
##       cat(paste("we appear to have a URL:",file,"\n"))
##     slashes <- gregexpr(pattern="/",text=file)[[1]]
##     if(length(slashes)<3)
##       cat(paste("dubious URL, it has only",length(slashes),"slashes\n"))
##     hostname <- substring(file,slashes[2]+1,slashes[3]-1)
##     if(debug)
##       cat(paste("hostname is",hostname,"\n"))

##     ## check that we can actually resolve the hostname
##     w <- options()$warn
##     options(warn=-1)
##     goodNet <- NULL
##     haveNSL <- exists("nsl")
##     if(haveNSL){
##       goodNet <- nsl(hostname)
##     else
      
##     if(debug)
##       if(is.null(goodNet))
##         cat(paste("nsl on",hostname,"returned NULL\n"))
##       else
##         cat(paste("nsl on",hostname,"returned",goodNet,"\n"))
##     if(is.null(goodNet)){
##       options(warn=w)
##       cat("Could not resolve the URL you provided.\n")
##       cat("Check the URL or your internet connection.\n")
##       return(invisible(NULL))
##     }
##     options(warn=w)
##   }

##   ## now actually try to read the data  
##   readResults <- try(readLines(file))
##   if("try-error" %in% class(readResults)){
##     cat(paste("readKH error: could not read from",file,"\n",
##               "execution terminating\n"))
##     data <- NULL
##   }
##   else{
##     data <- readResults
##     nRecs <- length(data)
##       cat(paste("read",file,"ok with",nRecs,"records\n"))
##     }
##   data
## }

dateExtract <- function(string){
  theMonths <- c("JANUARY","FEBRUARY","MARCH",
                 "APRIL","MAY","JUNE",
                 "JULY","AUGUST","SEPTEMBER",
                 "OCTOBER","NOVEMBER","DECEMBER")
  searchStringMonths <- paste(theMonths,collapse="|")
  foo <- unlist(strsplit(string,split=" "))
  foo <- foo[foo!=""]
  nFoo <- length(foo)
  whereMonth <- grep(pattern=searchStringMonths,foo)
  out <- ""
  if(length(whereMonth)==1){
    dateString <- foo[whereMonth:nFoo]
    dateString <- gsub(x=dateString,pattern=",",replacement="")
    ## dateString should be MONTH, DAY, YEAR
    month <- match(foo[whereMonth],table=theMonths)
    out <- paste(dateString[3],month,dateString[2],sep="-")
  }
  out
}

descriptionExtract <- function(recs){
  foo <- substring(recs,13)
  foo <- paste(foo,collapse="")
  foo <- gsub(foo,pattern="\n",replacement="")
  foo <- gsub(foo,pattern='[[:space:]]+',replacement=" ")
  foo <- gsub(foo,pattern='[[:space:]]$',replacement="")
  foo
}

## parse K&H dictionary files
dtlParser <- function(file,debug=TRUE){
  cat(paste("attempting to read dtl file",file,"\n"))
  warnLevel <- options()$warn
  options(warn=-1)
  data <- try(readLines(con=file),silent=TRUE)
  if(inherits(data,"try-error")){
    cat(paste("Could not read",file,"\n"))
    return(invisible(NULL))
  }
  options(warn=warnLevel)
  
  out <- NULL
  if(!is.null(data)){
    number <- as.numeric(substring(data,1,4))
    m <- unique(number)
    date <- rep(NA,length(m))
    description <- rep("",length(m))

    for(j in m){
      if(debug)
        cat(paste("dtlParser: processing record",j,"\n"))
      theRecs <- data[number==j]
      nRecs <- length(theRecs)

      ## extract the date
      date[j] <- dateExtract(theRecs[1])

      ## extract descriptive text
      if(nRecs>2)
        description[j] <- descriptionExtract(theRecs[3:nRecs])
    }
    out <- data.frame(date=date,
                      description=description)
    out$date <- as.Date(out$date,format="%Y-%m-%d")
  }

  out
}
