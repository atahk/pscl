\name{predict.zeroinfl}
\alias{predict.zeroinfl}
\alias{residuals.zeroinfl}
\alias{terms.zeroinfl}
\alias{model.matrix.zeroinfl}
\alias{coef.zeroinfl}
\alias{vcov.zeroinfl}
\alias{summary.zeroinfl}
\alias{print.summary.zeroinfl}
\alias{logLik.zeroinfl}
\alias{fitted.zeroinfl}
\alias{predprob.zeroinfl}
\alias{extractAIC.zeroinfl}

\title{Methods for zeroinfl Objects}
\description{
  Methods for extracting information from fitted zero-inflated
  regression model objects of class \code{"zeroinfl"}.
}

\usage{
\method{predict}{zeroinfl}(object, newdata,
  type = c("response", "prob", "count", "zero"), na.action = na.pass,
  at = NULL, \dots)
\method{residuals}{zeroinfl}(object, type = c("pearson", "response"), \dots)

\method{coef}{zeroinfl}(object, model = c("full", "count", "zero"), \dots)
\method{vcov}{zeroinfl}(object, model = c("full", "count", "zero"), \dots)

\method{terms}{zeroinfl}(x, model = c("count", "zero"), \dots)
\method{model.matrix}{zeroinfl}(object, model = c("count", "zero"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"zeroinfl"} as returned by
    \code{\link{zeroinfl}}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character specifying the type of predictions or residuals,
    respectively. For details see below.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{at}{optionally, if \code{type = "prob"}, a numeric vector at which
    the probabilities are evaluated. By default \code{0:max(y)} is used
    where \code{y} is the original observed response.}
  \item{model}{character specifying for which component of the model the
    terms or model matrix should be extracted.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"zeroinfl"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix
  (estimated from the Hessian in the numerical optimization of the log-likelihood).
  As usual, the \code{summary} method returns an object of class \code{"summary.zeroinfl"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.
  
  The methods for \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} by default
  return a single vector of coefficients and their associated covariance matrix,
  respectively, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model components can be extracted.
  
  Both the \code{\link[stats]{fitted}} and \code{\link[stats]{predict}} methods can
  compute fitted responses. The latter additionally provides the predicted density
  (i.e., probabilities for the observed counts), the predicted mean from the count
  component (without zero inflation) and the predicted probability for the zero
  component. The \code{\link[stats]{residuals}} method can compute
  raw residuals (observed - fitted) and Pearson residuals (raw residuals scaled by
  square root of variance function).
  
  The \code{\link[stats]{terms}} and \code{\link[stats]{model.matrix}} extractors can
  be used to extract the relevant information for either component of the model.

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{zeroinfl}}}

\examples{
data("bioChemists", package = "pscl")

fm_zip <- zeroinfl(art ~ ., data = bioChemists)
plot(residuals(fm_zip) ~ fitted(fm_zip))

coef(fm_zip)
coef(fm_zip, model = "count")

summary(fm_zip)
logLik(fm_zip)
}

\keyword{regression}
