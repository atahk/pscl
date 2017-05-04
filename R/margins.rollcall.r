## compute margins in a roll call object, add to the object in vote.data
computeMargins <- function(object,dropList=NULL){
  if(!inherits(object,"rollcall"))
    stop("margins only works on object of class rollcall.")

  tmpRollCall <- object
  if(!is.null(dropList)){
    tmpRollCall <- dropRollCall(object,dropList)
  }

  tab <- t(apply(convertCodes(tmpRollCall),
                 2,
                 marginfunc))
  tab <- cbind(tab,
               apply(tab,1,whichMinMargin))

  rownames(tab) <- dimnames(tmpRollCall$votes)[[2]]
  colnames(tab) <- c("Yea","Nay","NA","Min")
  tmpRollCall$voteMargins <- tab
    
  tmpRollCall
}

marginfunc <- function(x){
  ok <- !is.na(x)
  z <- c(sum(x[ok]==1),   ## Yeas
         sum(x[ok]==0),   ## Nays
         sum(!ok))        ## Missing
  z
}


marginWithCodes <- function(x,codes){
  n <- length(codes)
  tab <- rep(0,n)
  for(i in 1:n){
    if(is.list(codes))
      tab[i] <- sum(x %in% codes[[i]])
    else
      tab[i] <- sum(x %in% codes[i])
  }
  out <- c(tab,
           sum(tab),
           tab/sum(tab)*100)
  out
}


minMargin <- function(x){
  z <- rep(NA,2)
  ok <- !is.na(x)
  z[1] <- sum(x[ok]==0)
  z[2] <- sum(x[ok]==1)
  z <- z/sum(ok)
  out <- min(z)
  out
}

whichMinMargin <- function(x){
  x[which.min(x[1:2])]
}
