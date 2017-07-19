\name{predprob.ideal}
\alias{predprob.ideal}
\title{predicted probabilities from fitting ideal to rollcall data}
\description{
  Computes predicted probabilities of a \dQuote{Yea} vote conditional on
  the posterior means of the legislators' ideal points and
  vote-specific parameters.
}
\usage{
\method{predprob}{ideal}(obj, ...)
}
\arguments{
  \item{obj}{An object of class \code{\link{ideal}}}
  \item{\dots}{Arguments to be passed to other functions}
}
\details{
  This is a wrapper function to \code{\link{predict.ideal}}, extracting
  just the predicted probabilities component of the object returned by
  that function.  Predicted probabilities can and are generated for each
  voting decision, irrespective of whether the legislator actually voted
  on any particular roll call.
}
\value{
  A \code{\link{matrix}} of dimension \code{n} (number of legislators)
  by \code{m} (number of roll call votes).
}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\seealso{\code{\link{ideal}}, \code{\link{predprob}}, \code{\link{predict.ideal}}}
\examples{
f <- system.file("extdata","id1.rda",package="pscl")
load(f)
phat <- predprob(id1)
dim(phat)
}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{models}% __ONLY ONE__ keyword per line
