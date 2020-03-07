predprob <- function(obj, ...){
    UseMethod("predprob")
}

predprob.ideal <- function(obj, ...){
  if(!inherits(obj, "ideal"))
    stop("predprob.ideal only defined for objects of class ideal")
  else
    predict.ideal(obj,...)$pred.probs
}
