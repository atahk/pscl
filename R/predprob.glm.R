predprob.glm <- function(obj, newdata = NULL, at = NULL, ...){
    if(!inherits(obj,"glm"))
        stop("predprob.glm only available for glm objects (including class negbin)\n")

    isNegBin <- inherits(obj, "negbin")
    isPoisson <- family(obj)$family=="poisson"
    isBinomial <- family(obj)$family=="binomial"
    if(!isNegBin & !isPoisson & !isBinomial)
        stop(paste("your object of class",class(obj),"is unsupported by predprob.glm"))
    
    if(is.null(newdata))
        yhat <- predict(obj,
                        type="response")
    else
        yhat <- predict(obj,
                        newdata=newdata,
                        type="response")

    y <- obj$y
    yUnique <- if(is.null(at)) 0:max(y) else at
    nUnique <- length(yUnique)
    p <- matrix(NA,length(yhat),nUnique)
    dimnames(p) <- list(NULL,yUnique)

    if(isNegBin){
        for(i in 1:nUnique){
            p[,i] <- dnbinom(mu=yhat,
                             size=obj$theta,
                             x=yUnique[i])
        }
    }

    if(isPoisson){
        for(i in 1:nUnique){
            p[,i] <- dpois(lambda=yhat,
                           x=yUnique[i])
        }
    }

    if(isBinomial){
        if(is.null(newdata))
            p <- predict(obj,
                            type="response")
        else
            p <- predict(obj,
                            newdata=newdata,
                            type="response")
        p <- cbind(1-p,p)
        dimnames(p) <- list(NULL,c("0","1"))
    }

    p
    
}

