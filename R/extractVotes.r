## given an object of class ideal
## recover the rollcall object used in model fitting
## after applying the dropList etc

extractRollCallObject <- function(object){
  if(!inherits(object,"ideal"))
    stop("extractRollCallObject only defined for objects of class ideal")

  rcObj <- eval(object$call$object)
  dropList <- eval(object$call$dropList)

  tmpObj <- dropRollCall(rcObj,dropList)
  tmpObj
}

