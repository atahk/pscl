\name{predprob.glm}
\alias{predprob.glm}

\title{Predicted Probabilties for GLM Fits}

\description{
  Obtains predicted probabilities from a fitted generalized linear model
  object. 
}

\usage{
\method{predprob}{glm}(obj, newdata = NULL, at = NULL, ...)
}

\arguments{
  \item{obj}{a fitted object of class inheriting from \code{"glm"}}
  \item{newdata}{optionally, a data frame in which to look for variables
    with which to predict.  If omitted, the fitted linear predictors are
    used.}
  \item{at}{an optional numeric vector at which the probabilities are
    evaluated. By default \code{0:max(y)} where \code{y} is the original
    observed response.}
  \item{...}{arguments passed to or from other methods}
}

\details{
  This method is only defined for glm objects with
  \code{family=\link{binomial}} or \code{family=\link{poisson}}, or negative
  binomial count models fit with the \code{\link[MASS:glm.nb]{glm.nb}} function in \code{library(MASS)}.
}

\value{
  A matrix of predicted probabilities.  Each row in the matrix is a
  vector of probabilities, assigning predicted probabilities over the range of
  responses actually observed in the data.  For instance, for models
  with \code{family=binomial}, the matrix has two columns for the "zero"
  (or failure) and "one" (success) outcomes, respectively, and
  trivially, each row in the matrix sums to 1.0. 
  For counts fit with \code{family=poisson} or via \code{glm.nb}, the
  matrix has \code{length(0:max(y))} columns.  Each observation
  used in fitting the model generates a row to the returned matrix; alternatively, if
  \code{newdata} is supplied, the returned matrix will have as many rows
  as in \code{newdata}.
}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\seealso{\code{\link{predict.glm}}}
\examples{
data(bioChemists)
glm1 <- glm(art ~ .,
            data=bioChemists,
            family=poisson,
            trace=TRUE)  ## poisson GLM
phat <- predprob(glm1)
apply(phat,1,sum)                    ## almost all 1.0
}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{regression}% __ONLY ONE__ keyword per line
