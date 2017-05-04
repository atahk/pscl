\name{hurdletest}
\alias{hurdletest}

\title{Testing for the Presence of a Zero Hurdle}

\description{
  Wald test of the null hypothesis that no zero hurdle
  is required in hurdle regression models for count data.
}

\usage{
hurdletest(object, \dots)
}

\arguments{
  \item{object}{A fitted model object of class \code{"hurdle"}
    as returned by \code{\link{hurdle}}, see details for more information.}
  \item{\dots}{arguments passed to \code{\link[car]{linearHypothesis}}.}
}

\details{
  If the same count distribution and the same set of regressors is used in
  the hurdle model for both, the count component and the zero hurdle component,
  then a test of pairwise equality between all coefficients from the two
  components assesses the null hypothesis that no hurdle is needed in the
  model.
  
  The function \code{hurdletest} is a simple convenience interface to
  the function \code{\link[car]{linearHypothesis}} from the \pkg{car}
  packages that can be employed to carry out a Wald test for this 
  hypothesis.
}

\value{
  An object of class \code{"anova"} as returned by \code{\link[car]{linearHypothesis}}.
}

\references{
Cameron, A. Colin and Pravin K. Trivedi. 1998. \emph{Regression Analysis of Count 
Data}. New York: Cambridge University Press.

Cameron, A. Colin and Pravin K. Trivedi 2005. \emph{Microeconometrics: Methods and Applications}.
Cambridge: Cambridge University Press.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{hurdle}}, \code{\link[car]{linearHypothesis}}}

\examples{
data("bioChemists", package = "pscl")
fm <- hurdle(art ~ ., data = bioChemists, dist = "negbin", zero = "negbin")
hurdletest(fm)
}

\keyword{regression}
