\name{odTest}
\alias{odTest}

\title{likelihood ratio test for over-dispersion in count data}

\description{
  Compares the log-likelihoods of a negative binomial regression model
  and a Poisson regression model.
}

\usage{
odTest(glmobj, alpha=.05, digits = max(3, getOption("digits") - 3))
}

\arguments{
   \item{glmobj}{an object of class \code{negbin} produced by
     \code{\link[MASS:glm.nb]{glm.nb}}}  
   \item{alpha}{significance level of over-dispersion test}
   \item{digits}{ number of digits in printed output}
}

\details{The negative binomial model relaxes the assumption in the
  Poisson model that the (conditional) variance equals the (conditional)
  mean, by estimating one extra parameter.  A likelihood ratio (LR) test
  can be used to test the null hypothesis that the restriction implicit
  in the Poisson model is true. The LR test-statistic has a non-standard
  distribution, even asymptotically, since the negative binomial
  over-dispersion parameter (called \code{theta} in \code{glm.nb}) is
  restricted to be positive.  The asymptotic distribution of the LR
  (likelihood ratio) test-statistic has probability mass of one half at
  zero, and a half \eqn{\chi^2_1}{chi-square (1)} distribution above
  zero.  This means that if testing at the \eqn{\alpha}{alpha} = .05
  level, one should not reject the null unless the LR test statistic
  exceeds the critical value associated with the \eqn{2\alpha}{2 alpha}
  = .10 level; this LR test involves just one parameter restriction, so
  the critical value of the test statistic at the \eqn{p}{p} = .05 level
  is 2.7, instead of the usual 3.8 (i.e., the .90 quantile of the
  \eqn{\chi^2_1}{chi-square (1)} distribution, versus the .95 quantile).

  A Poisson model is run using \code{\link{glm}} with family set to
  \code{link{poisson}}, using the \code{\link{formula}} in the negbin
  model object passed as input.  The \code{\link{logLik}} functions are
  used to extract the log-likelihood for each model.
}

\value{None; prints results and returns silently}

\references{A. Colin Cameron and Pravin K. Trivedi (1998) \emph{Regression
    analysis of count data}. New York: Cambridge University Press.

  Lawless, J. F. (1987) Negative Binomial and Mixed Poisson
  Regressions. \emph{The Canadian Journal of Statistics}. 15:209-225.
 }

  \author{Simon Jackman \email{simon.jackman@sydney.edu.au}.  John Fox noted an
  error in an earlier version.}

\seealso{\code{\link[MASS:glm.nb]{glm.nb}}, \code{\link{logLik}}}

\examples{
data(bioChemists)
modelnb <- MASS::glm.nb(art ~ .,
                 data=bioChemists,
                 trace=TRUE)
odTest(modelnb)
}
\keyword{regression}% at least one, from doc/KEYWORDS
