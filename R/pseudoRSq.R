## pseudo r2 for binary, ordinal and multinomial models

pR2 <- function(object,...){
  UseMethod("pR2")
}

pR2Work <- function(llh,llhNull,n=NULL){
  if (is.null(n))
    n <- nobs(llh)
  McFadden <- 1 - llh/llhNull
  G2 <- -2*(llhNull-llh)
  r2ML <- 1 - exp(-G2/n)
  r2ML.max <- 1 - exp(llhNull*2/n)
  r2CU <- r2ML/r2ML.max
  out <- c(llh=llh,
           llhNull=llhNull,
           G2=G2,
           McFadden=McFadden,
           r2ML=r2ML,
           r2CU=r2CU)
  out      
}

pR2.default <- function(object,...){
  llh <- logLik(object)
  cat("fitting null model for pseudo-r2\n")
  objectNull <- update(object, ~ 1, data=model.frame(object))
  llhNull <- logLik(objectNull)
  pR2Work(llh,llhNull)
}
