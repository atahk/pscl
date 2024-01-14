\name{predict.hurdle}
\alias{predict.hurdle}
\alias{residuals.hurdle}
\alias{terms.hurdle}
\alias{model.matrix.hurdle}
\alias{coef.hurdle}
\alias{vcov.hurdle}
\alias{summary.hurdle}
\alias{print.summary.hurdle}
\alias{logLik.hurdle}
\alias{fitted.hurdle}
\alias{predprob.hurdle}
\alias{extractAIC.hurdle}

\title{Methods for hurdle Objects}
\description{
  Methods for extracting information from fitted hurdle
  regression model objects of class \code{"hurdle"}.
}

\usage{
\method{predict}{hurdle}(object, newdata,
  type = c("response", "prob", "count", "zero"), na.action = na.pass,
  at = NULL, \dots)
\method{residuals}{hurdle}(object, type = c("pearson", "response"), \dots)

\method{coef}{hurdle}(object, model = c("full", "count", "zero"), \dots)
\method{vcov}{hurdle}(object, model = c("full", "count", "zero"), \dots)

\method{terms}{hurdle}(x, model = c("count", "zero"), \dots)
\method{model.matrix}{hurdle}(object, model = c("count", "zero"), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"hurdle"} as returned by
    \code{\link{hurdle}}.}
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
  objects of class \code{"hurdle"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix
  (estimated from the Hessian in the numerical optimization of the log-likelihood).
  As usual, the \code{summary} method returns an object of class \code{"summary.hurdle"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.
  
  The methods for \code{\link[stats]{coef}} and \code{\link[stats]{vcov}} by default
  return a single vector of coefficients and their associated covariance matrix,
  respectively, i.e., all coefficients are concatenated. By setting the \code{model}
  argument, the estimates for the corresponding model component can be extracted.
  
  Both the \code{\link[stats]{fitted}} and \code{\link[stats]{predict}} methods can
  compute fitted responses. The latter additionally provides the predicted density
  (i.e., probabilities for the observed counts), the predicted mean from the count
  component (without zero hurdle) and the predicted ratio of probabilities for
  observing a non-zero count. The latter is the ratio of probabilities for a non-zero
  implied by the zero hurdle component and a non-zero count in the non-truncated
  count distribution. See also Appendix C in Zeileis et al. (2008).
  
  The \code{\link[stats]{residuals}} method can compute raw residuals
  (observed - fitted) and Pearson residuals (raw residuals scaled by
  square root of variance function).
  
  The \code{\link[stats]{terms}} and \code{\link[stats]{model.matrix}} extractors can
  be used to extract the relevant information for either component of the model.

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\references{
  Zeileis, Achim, Christian Kleiber and Simon Jackman 2008.
  \dQuote{Regression Models for Count Data in R.} 
  \emph{Journal of Statistical Software}, \bold{27}(8).
  URL \url{https://www.jstatsoft.org/v27/i08/}.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{hurdle}}}

\examples{
data("bioChemists", package = "pscl")
fm <- hurdle(art ~ ., data = bioChemists)

plot(residuals(fm) ~ fitted(fm))

coef(fm)
coef(fm, model = "zero")

summary(fm)
logLik(fm)
}

\keyword{regression}
