predprob <- function(obj, ...){
    UseMethod("predprob")
}

predprob.ideal <- function(obj, ...){
  if(!("ideal" %in% class(obj)))
    stop("predprob.ideal only defined for objects of class ideal")
  else
    predict.ideal(obj,...)$pred.probs
}
