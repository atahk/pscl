hurdle <- function(formula, data, subset, na.action, weights, offset,
                   dist = c("poisson", "negbin", "geometric"),
                   zero.dist = c("binomial", "poisson", "negbin", "geometric"),
                   link = c("logit", "probit", "cloglog", "cauchit", "log"),
		   control = hurdle.control(...),
		   model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up likelihood components
  zeroPoisson <- function(parms) {
    ## mean
    mu <- as.vector(exp(Z %*% parms + offsetz))
    ## log-likelihood
    loglik0 <- -mu ## = dpois(0, lambda = mu, log = TRUE)
    ## collect and return
    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * log(1 - exp(loglik0[Y1])))
    loglik
  }

  countPoisson <- function(parms) {
    ## mean
    mu <- as.vector(exp(X %*% parms + offsetx))[Y1]
    ## log-likelihood
    loglik0 <- -mu ## = dpois(0, lambda = mu, log = TRUE)
    loglik1 <- dpois(Y[Y1], lambda = mu, log = TRUE)
    ## collect and return
    loglik <- sum(weights[Y1] * loglik1) - sum(weights[Y1] * log(1 - exp(loglik0)))
    loglik
  }

  zeroNegBin <- function(parms) {
    ## parameters
    mu <- as.vector(exp(Z %*% parms[1:kz] + offsetz))
    theta <- exp(parms[kz+1])
    ## log-likelihood
    loglik0 <- suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))
    ## collect and return
    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * log(1 - exp(loglik0[Y1])))
    loglik
  }

  countNegBin <- function(parms) {
    ## parameters
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))[Y1]
    theta <- exp(parms[kx+1])
    ## log-likelihood
    loglik0 <- suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))
    loglik1 <- suppressWarnings(dnbinom(Y[Y1], size = theta, mu = mu, log = TRUE))
    ## collect and return
    loglik <- sum(weights[Y1] * loglik1) - sum(weights[Y1] * log(1 - exp(loglik0)))
    loglik
  }

  zeroGeom <- function(parms) zeroNegBin(c(parms, 0))
  
  countGeom <- function(parms) countNegBin(c(parms, 0))

  zeroBinom <- function(parms) {
    ## mean
    mu <- as.vector(linkinv(Z %*% parms + offsetz))
    ## log-likelihood
    loglik <- sum(weights[Y0] * log(1 - mu[Y0])) + sum(weights[Y1] * log(mu[Y1]))
    loglik  
  }

  countGradPoisson <- function(parms) {
    eta <- as.vector(X %*% parms + offsetx)[Y1]
    mu <- exp(eta)
    colSums(((Y[Y1] - mu) - exp(ppois(0, lambda = mu, log.p = TRUE) -
      ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)) * weights[Y1] * X[Y1, , drop = FALSE])
  }
  
  countGradGeom <- function(parms) {
    eta <- as.vector(X %*% parms + offsetx)[Y1]
    mu <- exp(eta)      
    colSums(((Y[Y1] - mu * (Y[Y1] + 1)/(mu + 1)) -
      exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
        pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) -
	log(mu + 1) + eta)) * weights[Y1] * X[Y1, , drop = FALSE])
  }

  countGradNegBin <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)[Y1]
    mu <- exp(eta)      
    theta <- exp(parms[kx+1])
    logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
        pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
    rval <- colSums(((Y[Y1] - mu * (Y[Y1] + theta)/(mu + theta)) -
      exp(logratio + log(theta) - log(mu + theta) + eta)) * weights[Y1] * X[Y1, , drop = FALSE])
    rval2 <- sum((digamma(Y[Y1] + theta) - digamma(theta) +    
      log(theta) - log(mu + theta) + 1 - (Y[Y1] + theta)/(mu + theta) +
      exp(logratio) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))) * weights[Y1]) * theta
    c(rval, rval2)
  }  
  
  zeroGradPoisson <- function(parms) {
    eta <- as.vector(Z %*% parms + offsetz)
    mu <- exp(eta)
    colSums(ifelse(Y0, -mu, exp(ppois(0, lambda = mu, log.p = TRUE) -
      ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE) + eta)) * weights * Z)
  }

  zeroGradGeom <- function(parms) {
    eta <- as.vector(Z %*% parms + offsetz)
    mu <- exp(eta)
    colSums(ifelse(Y0, -mu/(mu + 1), exp(pnbinom(0, mu = mu, size = 1, log.p = TRUE) -
      pnbinom(0, mu = mu, size = 1, lower.tail = FALSE, log.p = TRUE) - log(mu + 1) + eta)) * weights * Z)
  }

  zeroGradNegBin <- function(parms) {
    eta <- as.vector(Z %*% parms[1:kz] + offsetz)
    mu <- exp(eta)
    theta <- exp(parms[kz+1])
    logratio <- pnbinom(0, mu = mu, size = theta, log.p = TRUE) -
      pnbinom(0, mu = mu, size = theta, lower.tail = FALSE, log.p = TRUE)
    rval <- colSums(ifelse(Y0, -mu * theta/(mu + theta),
      exp(logratio + log(theta) - log(mu + theta) + eta)) * weights * Z)
    rval2 <- sum(ifelse(Y0, log(theta) - log(mu + theta) + 1 - theta/(mu + theta),
      -exp(logratio) * (log(theta) - log(mu + theta) + 1 - theta/(mu + theta))) * weights * theta)
    c(rval, rval2)
  }

  zeroGradBinom <- function(parms) {
    eta <- as.vector(Z %*% parms + offsetz)
    mu <- linkinv(eta)
    colSums(ifelse(Y0, -1/(1-mu), 1/mu)  * linkobj$mu.eta(eta) * weights * Z)  
  }

  ## collect likelihood components
  dist <- match.arg(dist)
  zero.dist <- match.arg(zero.dist)
  countDist <- switch(dist,
                     "poisson" = countPoisson,
		     "geometric" = countGeom,
		     "negbin" = countNegBin)
  zeroDist <- switch(zero.dist,
                     "poisson" = zeroPoisson,
		     "geometric" = zeroGeom,
		     "negbin" = zeroNegBin,
		     "binomial" = zeroBinom)
  countGrad <- switch(dist,
                     "poisson" = countGradPoisson,
		     "geometric" = countGradGeom,
		     "negbin" = countGradNegBin)
  zeroGrad <- switch(zero.dist,
                     "poisson" = zeroGradPoisson,
		     "geometric" = zeroGradGeom,
		     "negbin" = zeroGradNegBin,
		     "binomial" = zeroGradBinom)
  loglikfun <- function(parms) countDist(parms[1:(kx + (dist == "negbin"))]) +
    zeroDist(parms[(kx + (dist == "negbin") + 1):(kx + kz + (dist == "negbin") + (zero.dist == "negbin"))])
  gradfun <- function(parms) c(countGrad(parms[1:(kx + (dist == "negbin"))]),
    zeroGrad(parms[(kx + (dist == "negbin") + 1):(kx + kz + (dist == "negbin") + (zero.dist == "negbin"))]))


  ## binary link processing
  linkstr <- match.arg(link)
  linkobj <- make.link(linkstr)
  linkinv <- linkobj$linkinv

  if(control$trace) cat("Hurdle Count Model\n",
    paste("count model:", dist, "with log link\n"),
    paste("zero hurdle model:", zero.dist, "with", ifelse(zero.dist == "binomial", linkstr, "log"), "link\n"),
    sep = "")


  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## extended formula processing
  if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
  {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~ .
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  } else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if(inherits(try(terms(ffz), silent = TRUE), "try-error")) {
    ffz <- eval(parse(text = sprintf( paste("%s -", deparse(ffc[[2]])), deparse(ffz) )))
  }

  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrices, response
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~ .), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")


  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(all(Y > 0)) stop("invalid dependent variable, minimum count is not zero")  
  if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
    stop("invalid dependent variable, non-integer values")
  Y <- as.integer(round(Y + 0.001))
  if(any(Y < 0)) stop("invalid dependent variable, negative counts")
  if(zero.dist == "negbin" & isTRUE(all.equal(as.vector(Z), rep.int(Z[1], length(Z)))))
    stop("negative binomial zero hurdle model is not identified with only an intercept")
  
  if(control$trace) {
    cat("dependent variable:\n")
    tab <- table(factor(Y, levels = 0:max(Y)), exclude = NULL)
    names(dimnames(tab)) <- NULL
    print(tab)
  }

  ## convenience variables
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- Y <= 0
  Y1 <- Y > 0

  ## weights and offset
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
  if(is.null(offsetx)) offsetx <- 0
  if(length(offsetx) == 1) offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
  if(is.null(offsetz)) offsetz <- 0
  if(length(offsetz) == 1) offsetz <- rep.int(offsetz, n)
  offsetz <- as.vector(offsetz)

  ## starting values
  start <- control$start
  if(!is.null(start)) {
    valid <- TRUE
    if(!("count" %in% names(start))) {
      valid <- FALSE
      warning("invalid starting values, count model coefficients not specified")
      start$count <- rep.int(0, kx)
    }
    if(!("zero" %in% names(start))) {
      valid <- FALSE
      warning("invalid starting values, zero-inflation model coefficients not specified")
      start$zero <- rep.int(0, kz)
    }
    if(length(start$count) != kx) {
      valid <- FALSE
      warning("invalid starting values, wrong number of count model coefficients")
    }
    if(length(start$zero) != kz) {
      valid <- FALSE
      warning("invalid starting values, wrong number of zero-inflation model coefficients")
    }
    if(dist == "negbin" | zero.dist == "negbin") {
      if(!("theta" %in% names(start))) start$theta <- c(1, 1)
      start <- list(count = start$count, zero = start$zero, theta = rep(start$theta, length.out = 2))
      if(is.null(names(start$theta))) names(start$theta) <- c("count", "zero")      
      if(dist != "negbin") start$theta <- start$theta["zero"]
      if(zero.dist != "negbin") start$theta <- start$theta["count"]
    } else {
      start <- list(count = start$count, zero = start$zero)
    }
    if(!valid) start <- NULL
  }
  
  if(is.null(start)) {
    if(control$trace) cat("generating starting values...")
    model_count <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)
    model_zero <- switch(zero.dist,
      "poisson" = glm.fit(Z, Y, family = poisson(), weights = weights, offset = offsetz),
      "negbin" = glm.fit(Z, Y, family = poisson(), weights = weights, offset = offsetz),
      "geometric" = suppressWarnings(glm.fit(Z, factor(Y > 0), family = binomial(), weights = weights, offset = offsetz)),
      "binomial" = suppressWarnings(glm.fit(Z, factor(Y > 0), family = binomial(link = linkstr), weights = weights, offset = offsetz)))
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
    start$theta <- c(count = if(dist == "negbin") 1 else NULL,
                     zero = if(zero.dist == "negbin") 1 else NULL)
    if(control$trace) cat("done\n")
  }

  ## model fitting
  ## control parameters
  method <- control$method
  hessian <- control$hessian
  separate <- control$separate
  ##ocontrol <- control
  control$method <- control$hessian <- control$separate <- control$start <- NULL

  ## ML estimation
  ## separate estimation of censored and truncated component...
  if(separate) {
    if(control$trace) cat("calling optim() for count component estimation:\n")
    fit_count <- optim(fn = countDist, gr = countGrad,
      par = c(start$count, if(dist == "negbin") log(start$theta["count"]) else NULL),
      method = method, hessian = hessian, control = control)

    if(control$trace) cat("calling optim() for zero hurdle component estimation:\n")
    fit_zero <- optim(fn = zeroDist, gr = zeroGrad,
      par = c(start$zero,  if(zero.dist == "negbin") log(start$theta["zero"]) else NULL),
      method = method, hessian = hessian, control = control)
    if(control$trace) cat("done\n")
    fit <- list(count = fit_count, zero = fit_zero)

    ## coefficients
    coefc <- fit_count$par[1:kx]
    coefz <- fit_zero$par[1:kz]
    theta <- c(count = if(dist == "negbin") as.vector(exp(fit_count$par[kx+1])) else NULL,
               zero = if(zero.dist == "negbin") as.vector(exp(fit_zero$par[kz+1])) else NULL)
    ## covariances
    vc_count <- tryCatch(-solve(as.matrix(fit_count$hessian)),
                         error=function(e) {
                             warning(e$message, call=FALSE)
                             k <- nrow(as.matrix(fit_count$hessian))
                             return(matrix(NA, k, k))
                         })
    vc_zero <- tryCatch(-solve(as.matrix(fit_zero$hessian)),
                        error=function(e) {
                            warning(e$message, call=FALSE)
                            k <- nrow(as.matrix(fit_zero$hessian))
                            return(matrix(NA, k, k))
                        })
    SE.logtheta <- list()
    if(dist == "negbin") {
      SE.logtheta$count <- as.vector(sqrt(diag(vc_count)[kx+1]))
      vc_count <- vc_count[-(kx+1), -(kx+1), drop = FALSE]
    }
    if(zero.dist == "negbin") {
      SE.logtheta$zero <- as.vector(sqrt(diag(vc_zero)[kz+1]))
      vc_zero <- vc_zero[-(kz+1), -(kz+1), drop = FALSE]
    }
    vc <- rbind(cbind(vc_count, matrix(0, kx, kz)), cbind(matrix(0, kz, kx), vc_zero))
    SE.logtheta <- unlist(SE.logtheta)
  } else {
  ## ...or joint.
    if(control$trace) cat("calling optim() for joint count and zero hurlde estimation:\n")
    fit <- optim(fn = loglikfun, gr = gradfun,
      par = c(start$count, if(dist == "negbin") log(start$theta["count"]) else NULL,
              start$zero,  if(zero.dist == "negbin") log(start$theta["zero"]) else NULL),
      method = method, hessian = hessian, control = control)
    if(fit$convergence > 0) warning("optimization failed to converge")
    if(control$trace) cat("done\n")

    ## coefficients
    coefc <- fit$par[1:kx]
    coefz <- fit$par[(kx + (dist == "negbin") + 1):(kx + kz + (dist == "negbin"))]
    ## covariances
    vc <- tryCatch(-solve(as.matrix(fit$hessian)),
                   error=function(e) {
                       warning(e$message, call=FALSE)
                       k <- nrow(as.matrix(fit$hessian))
                       return(matrix(NA, k, k))
                   })
    np <- c(if(dist == "negbin") kx+1 else NULL,
            if(zero.dist == "negbin") kx+kz+1+(dist == "negbin") else NULL)
    if(length(np) > 0) {
      theta <- as.vector(exp(fit$par[np]))
      SE.logtheta <- as.vector(sqrt(diag(vc)[np]))
      names(theta) <- names(SE.logtheta) <- c(if(dist == "negbin") "count" else NULL,
                                              if(zero.dist == "negbin") "zero" else NULL)
      vc <- vc[-np, -np, drop = FALSE]
    } else {
      theta <- NULL
      SE.logtheta <- NULL
    }

  }
  names(coefc) <- names(start$count) <- colnames(X)
  names(coefz) <- names(start$zero) <- colnames(Z)
  colnames(vc) <- rownames(vc) <- c(paste("count", colnames(X), sep = "_"),
                                    paste("zero",  colnames(Z), sep = "_"))

  ## fitted and residuals
  phi <- if(zero.dist == "binomial") linkinv(Z %*% coefz + offsetz)[,1] else exp(Z %*% coefz + offsetz)[,1]
  p0_zero <- switch(zero.dist,
  		    "binomial" = log(phi),
        	    "poisson" = ppois(0, lambda = phi, lower.tail = FALSE, log.p = TRUE),
        	    "negbin" = pnbinom(0, size = theta["zero"], mu = phi, lower.tail = FALSE, log.p = TRUE),
        	    "geometric" = pnbinom(0, size = 1, mu = phi, lower.tail = FALSE, log.p = TRUE))

  mu <- exp(X %*% coefc + offsetx)[,1]
  p0_count <- switch(dist,
        	    "poisson" = ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE),
        	    "negbin" = pnbinom(0, size = theta["count"], mu = mu, lower.tail = FALSE, log.p = TRUE),
        	    "geometric" = pnbinom(0, size = 1, mu = mu, lower.tail = FALSE, log.p = TRUE))
  Yhat <- exp((p0_zero - p0_count) + log(mu))

  res <- sqrt(weights) * (Y - Yhat)

  ## effective observations
  nobs <- sum(weights > 0) ## = n - sum(weights == 0)

  rval <- list(coefficients = list(count = coefc, zero = coefz),
    residuals = res,
    fitted.values = Yhat,
    optim = fit,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep.int(1L, n))) NULL else weights,
    offset = list(count = if(identical(offsetx, rep.int(0, n))) NULL else offsetx,
      zero = if(identical(offsetz, rep.int(0, n))) NULL else offsetz),
    n = nobs,
    df.null = nobs - 2,
    df.residual = nobs - (kx + kz + (dist == "negbin") + (zero.dist == "negbin")),
    terms = list(count = mtX, zero = mtZ, full = mt),
    theta = theta,
    SE.logtheta = SE.logtheta,
    loglik = if(separate) fit_count$value + fit_zero$value else fit$value,
    vcov = vc,
    dist = list(count = dist, zero = zero.dist),
    link = if(zero.dist == "binomial") linkstr else NULL,
    linkinv = if(zero.dist == "binomial") linkinv else NULL,
    separate = separate,
    converged = if(separate) fit_count$convergence < 1 & fit_zero$convergence < 1 else fit$convergence < 1,
    call = cl,
    formula = ff,
    levels = .getXlevels(mt, mf),
    contrasts = list(count = attr(X, "contrasts"), zero = attr(Z, "contrasts"))
  )
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(count = X, zero = Z)
      
  class(rval) <- "hurdle"
  return(rval)
}

hurdle.control <- function(method = "BFGS", maxit = 10000, trace = FALSE, separate = TRUE, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, trace = trace, separate = separate, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(!is.null(rval$hessian)) warning("hessian must not be modified")
  rval$hessian <- TRUE
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.6)
  rval
}

coef.hurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  rval <- object$coefficients
  rval <- switch(model,
                 "full" = structure(c(rval$count, rval$zero),
		   .Names = c(paste("count", names(rval$count), sep = "_"),
                   paste("zero", names(rval$zero), sep = "_"))),
		 "count" = rval$count,
		 "zero" = rval$zero)
  rval
}

vcov.hurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  rval <- object$vcov
  if(model == "full") return(rval)

  cf <- object$coefficients[[model]]
  wi <- seq(along = object$coefficients$count)
  rval <- if(model == "count") rval[wi, wi, drop = FALSE] else rval[-wi, -wi, drop = FALSE]
  colnames(rval) <- rownames(rval) <- names(cf)
  return(rval)
}

logLik.hurdle <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, nobs = object$n, class = "logLik")
}

print.hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Count model coefficients (truncated ", x$dist$count, " with log link):\n", sep = ""))
    print.default(format(x$coefficients$count, digits = digits), print.gap = 2, quote = FALSE)
    if(x$dist$count == "negbin") cat(paste("Theta =", round(x$theta["count"], digits), "\n"))
  
    zero_dist <- if(x$dist$zero != "binomial") paste("censored", x$dist$zero, "with log link")
      else paste("binomial with", x$link, "link")
    cat(paste("\nZero hurdle model coefficients (", zero_dist, "):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits), print.gap = 2, quote = FALSE)
    if(x$dist$zero == "negbin") cat(paste("Theta =", round(x$theta["zero"], digits), "\n"))
    cat("\n")
  }
  
  invisible(x)
}

summary.hurdle <- function(object,...)
{
  ## residuals
  object$residuals <- residuals(object, type = "pearson")
  
  ## compute z statistics
  kc <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients$count, object$coefficients$zero)  
  if(object$dist$count == "negbin") {
    coef <- c(coef[1:kc], "Log(theta)" = as.vector(log(object$theta["count"])), coef[(kc+1):(kc+kz)])
    se <- c(se[1:kc], object$SE.logtheta["count"], se[(kc+1):(kc+kz)])
    kc <- kc+1
  }
  if(object$dist$zero == "negbin") {
    coef <- c(coef, "Log(theta)" = as.vector(log(object$theta["zero"])))
    se <- c(se, object$SE.logtheta["zero"])
    kz <- kz+1
  }
  zstat <- coef/se
  pval <- 2*pnorm(-abs(zstat))
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients$count <- coef[1:kc,,drop = FALSE]
  object$coefficients$zero <- coef[(kc+1):(kc+kz),,drop = FALSE]
  
  ## number of iterations
  object$iterations <- if(!object$separate) tail(na.omit(object$optim$count), 1)
    else tail(na.omit(object$optim$count$count), 1) + tail(na.omit(object$optim$zero$count), 1)
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- object$separate <- NULL

  ## return
  class(object) <- "summary.hurdle"
  object
}

print.summary.hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {

    cat("Pearson residuals:\n")
    print(structure(quantile(x$residuals),
      names = c("Min", "1Q", "Median", "3Q", "Max")), digits = digits, ...)  

    cat(paste("\nCount model coefficients (truncated ", x$dist$count, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
  
    zero_dist <- if(x$dist$zero != "binomial") paste("censored", x$dist$zero, "with log link")
      else paste("binomial with", x$link, "link")
    cat(paste("Zero hurdle model coefficients (", zero_dist, "):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    
    if(getOption("show.signif.stars") & any(rbind(x$coefficients$count, x$coefficients$zero)[,4] < 0.1, na.rm=TRUE))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    if(!is.null(x$theta)) cat(paste("\nTheta:", paste(names(x$theta), round(x$theta, digits), sep = " = ", collapse = ", ")))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }
  
  invisible(x)
}

terms.hurdle <- function(x, model = c("count", "zero"), ...) {
  x$terms[[match.arg(model)]]
}

model.matrix.hurdle <- function(object, model = c("count", "zero"), ...) {
  model <- match.arg(model)
  if(!is.null(object$x)) rval <- object$x[[model]]
    else if(!is.null(object$model)) rval <- model.matrix(object$terms[[model]], object$model, contrasts = object$contrasts[[model]])
    else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

predict.hurdle <- function(object, newdata, type = c("response", "prob", "count", "zero"),
  na.action = na.pass, at = NULL, ...)
{
    type <- match.arg(type)

    ## if no new data supplied
    if(missing(newdata)) {
      if(type != "response") {
        if(!is.null(object$x)) {
	  X <- object$x$count
	  Z <- object$x$zero
	} else if(!is.null(object$model)) {
          X <- model.matrix(object$terms$count, object$model, contrasts = object$contrasts$count)
          Z <- model.matrix(object$terms$zero,  object$model, contrasts = object$contrasts$zero)	
	} else {
	  stop("predicted probabilities cannot be computed with missing newdata")
	}
	offsetx <- if(is.null(object$offset$count)) rep.int(0, NROW(X)) else object$offset$count
	offsetz <- if(is.null(object$offset$zero))  rep.int(0, NROW(Z)) else object$offset$zero
      } else {
        return(object$fitted.values)
      }
    } else {
      mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels)
      X <- model.matrix(delete.response(object$terms$count), mf, contrasts = object$contrasts$count)
      Z <- model.matrix(delete.response(object$terms$zero),  mf, contrasts = object$contrasts$zero)
      offsetx <- model_offset_2(mf, terms = object$terms$count, offset = FALSE)
      offsetz <- model_offset_2(mf, terms = object$terms$zero,  offset = FALSE)
      if(is.null(offsetx)) offsetx <- rep.int(0, NROW(X))
      if(is.null(offsetz)) offsetz <- rep.int(0, NROW(Z))
      if(!is.null(object$call$offset)) offsetx <- offsetx + eval(object$call$offset, newdata)
    }

    phi <- if(object$dist$zero == "binomial") object$linkinv(Z %*% object$coefficients$zero + offsetz)[,1]
      else exp(Z %*% object$coefficients$zero + offsetz)[,1]
    p0_zero <- switch(object$dist$zero,
                      "binomial" = log(phi),
		      "poisson" = ppois(0, lambda = phi, lower.tail = FALSE, log.p = TRUE),
		      "negbin" = pnbinom(0, size = object$theta["zero"], mu = phi, lower.tail = FALSE, log.p = TRUE),
		      "geometric" = pnbinom(0, size = 1, mu = phi, lower.tail = FALSE, log.p = TRUE))

    mu <- exp(X %*% object$coefficients$count + offsetx)[,1]
    p0_count <- switch(object$dist$count,
		      "poisson" = ppois(0, lambda = mu, lower.tail = FALSE, log.p = TRUE),
		      "negbin" = pnbinom(0, size = object$theta["count"], mu = mu, lower.tail = FALSE, log.p = TRUE),
		      "geometric" = pnbinom(0, size = 1, mu = mu, lower.tail = FALSE, log.p = TRUE))
    logphi <- p0_zero - p0_count
    
    if(type == "response") rval <- exp(logphi + log(mu))
    if(type == "count") rval <- mu
    if(type == "zero") rval <- exp(logphi)

    ## predicted probabilities
    if(type == "prob") {
      if(!is.null(object$y)) y <- object$y
        else if(!is.null(object$model)) y <- model.response(object$model)
	else stop("predicted probabilities cannot be computed for fits with y = FALSE and model = FALSE")

      yUnique <- if(is.null(at)) 0:max(y) else at
      nUnique <- length(yUnique)
      rval <- matrix(NA, nrow = length(mu), ncol = nUnique)
      dimnames(rval) <- list(rownames(X), yUnique)
      
      rval[,1] <- 1 - exp(p0_zero)
      switch(object$dist$count,
             "poisson" = {
               for(i in 2:nUnique) rval[,i] <- exp(logphi + dpois(yUnique[i], lambda = mu, log = TRUE))
	     },
	     "negbin" = {
               for(i in 2:nUnique) rval[,i] <- exp(logphi + dnbinom(yUnique[i], mu = mu, size = object$theta["count"], log = TRUE))
	     },
	     "geometric" = {
               for(i in 2:nUnique) rval[,i] <- exp(logphi + dnbinom(yUnique[i], mu = mu, size = 1, log = TRUE))
	     })
    }
   
    rval
}

fitted.hurdle <- function(object, ...) {
  object$fitted.values
}

residuals.hurdle <- function(object, type = c("pearson", "response"), ...) {

  type <- match.arg(type)
  res <- object$residuals

  switch(type,
  
  "response" = {
    return(res)
  },
  
  "pearson" = {
    mu <- predict(object, type = "count")
    phi <- predict(object, type = "zero")
    theta1 <- switch(object$dist$count,
      "poisson" = 0,
      "geometric" = 1,
      "negbin" = 1/object$theta["count"])
    vv <- object$fitted.values * (1 + ((1-phi) + theta1) * mu)
    return(res/sqrt(vv))
  })
}

predprob.hurdle <- function(obj, ...){
    predict(obj, type = "prob", ...)
}

extractAIC.hurdle <- function(fit, scale = NULL, k = 2, ...) {
  c(attr(logLik(fit), "df"), AIC(fit, k = k))
}

hurdletest <- function(object, ...) {
  stopifnot(inherits(object, "hurdle"))
  stopifnot(object$dist$count == object$dist$zero)
  stopifnot(all(sort(names(object$coefficients$count)) == sort(names(object$coefficients$zero))))
  stopifnot(requireNamespace("car"))
  nam <- names(object$coefficients$count)
  lh <- paste("count_", nam, " = ", "zero_", nam, sep = "")
  rval <- car::linearHypothesis(object, lh, ...)
  attr(rval, "heading")[1] <- "Wald test for hurdle models\n\nRestrictions:"
  return(rval)
}

## convenience helper function
model_offset_2 <- function(x, terms = NULL, offset = TRUE)
## allow optionally different terms
## potentially exclude "(offset)"
{
  if(is.null(terms)) terms <- attr(x, "terms")
  offsets <- attr(terms, "offset")
  if(length(offsets) > 0) {
    ans <- if(offset) x$"(offset)" else NULL
    if(is.null(ans)) ans <- 0
    for(i in offsets) ans <- ans + x[[deparse(attr(terms, "variables")[[i + 1]])]]
    ans
  }
  else {
    ans <- if(offset) x$"(offset)" else NULL
  }
  if(!is.null(ans) && !is.numeric(ans)) stop("'offset' must be numeric")
  ans
}
