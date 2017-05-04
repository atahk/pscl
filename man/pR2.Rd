\name{pR2}
\alias{pR2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compute various pseudo-R2 measures}
\description{
  compute various pseudo-R2 measures for various GLMs
}
\usage{
pR2(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model object, for now of class \code{glm},
    \code{polr}, or \code{mulitnom}} 
  \item{\dots}{additional arguments to be passed to or from functions}
}

\details{Numerous pseudo r-squared measures have been proposed for
  generalized linear models, involving a comparison of the
  log-likelihood for the fitted model against the log-likelihood of a
  null/restricted model with no predictors, normalized to run from zero
  to one as the fitted model provides a better fit to the data
  (providing a rough analogue to the computation of r-squared in a
  linear regression).  }

\value{
  A vector of length 6 containing
  \item{llh}{The log-likelihood from the fitted model}
  \item{llhNull}{The log-likelihood from the intercept-only restricted
    model}
  \item{G2}{Minus two times the difference in the log-likelihoods}
  \item{McFadden}{McFadden's pseudo r-squared}
  \item{r2ML}{Maximum likelihood pseudo r-squared}
  \item{r2CU}{Cragg and Uhler's pseudo r-squared}
}
\references{Long, J. Scott. 1997. \emph{Regression Models for Categorical and
  Limited Dependent Variables}.  Sage.  pp104-106.}
\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{extractAIC}}, \code{\link{logLik}}}
\examples{
data(admit)
## ordered probit model
op1 <- MASS::polr(score ~ gre.quant + gre.verbal + ap + pt + female,
            Hess=TRUE,
            data=admit,
            method="probit")
pR2(op1)   
}
\keyword{models}
