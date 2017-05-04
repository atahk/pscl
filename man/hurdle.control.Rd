\name{hurdle.control}
\alias{hurdle.control}

\title{Control Parameters for Hurdle Count Data Regression}
\description{
  Various parameters that control fitting of hurdle regression models
  using \code{\link{hurdle}}.
}

\usage{
hurdle.control(method = "BFGS", maxit = 10000, trace = FALSE,
  separate = TRUE, start = NULL, \dots)
}

\arguments{
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{trace}{logical or integer controlling whether tracing information on  
    the progress of the optimization should be produced (passed to \code{\link[stats]{optim}}).}
  \item{separate}{logical. Should the estimation of the parameters in the truncated
    count component and hurdle zero component be carried out separately? See details.}
  \item{start}{an optional list with elements \code{"count"} and \code{"zero"}
    (and potentially \code{"theta"}) containing the coefficients for the corresponding component.}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{hurdle}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{hurdle.control}}.
  Most arguments are passed on directly to \code{optim}, only \code{trace} is also
  used within \code{hurdle} and \code{separate}/\code{start} control how \code{optim}
  is called.

  Starting values can be supplied via \code{start} or estimated by
  \code{\link[stats]{glm.fit}} (default).
  If \code{separate = TRUE} (default) the likelihoods of the truncated count component and the
  hurdle zero component will be maximized separately, otherwise the joint likelihood is
  set up and maximized. Standard errors are derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. To supply starting 
  values, \code{start} should be a list with elements \code{"count"} and \code{"zero"}
  and potentially \code{"theta"} (a named vector, for models with negative binomial components 
  only) containing the starting values for the coefficients of the corresponding component 
  of the model.
}

\value{
  A list with the arguments specified.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{hurdle}}}

\examples{
data("bioChemists", package = "pscl")

## default start values
fm1 <- hurdle(art ~ fem + ment, data = bioChemists,
              dist = "negbin", zero = "negbin")

## user-supplied start values and other options
fm2 <- hurdle(art ~ fem + ment, data = bioChemists,
              dist = "negbin",
              zero = "negbin",
              trace=TRUE,
              separate=FALSE,
              start = list(count = c(0.3, -0.2, 0),
                           zero = c(4, -2, 0.8),
                           theta = c(count = 2, zero = 0.1)))

}

\keyword{regression}
