\name{plot.predict.ideal}
\alias{plot.predict.ideal}

\title{plot methods for predictions from ideal objects}

\description{
  Plot classification success rates by legislators, or by roll calls,
  using predictions from ideal.
}

\usage{
\method{plot}{predict.ideal}(x, type = c("legis", "votes"),...)
}

\arguments{
  \item{x}{an object of class \code{\link{predict.ideal}}.}
  \item{type}{string; one of \code{legis} or \code{votes}.}
  \item{...}{further arguments passed to or from other methods.}
}

\details{
  \code{type="legis"} produces a plot of the \dQuote{percent correctly
    predicted} for each legislator/subject
  (using the classification threshold set in
  \code{\link{predict.ideal}}) against the estimated ideal point of each
  legislator/subject (the estimated mean of the posterior density of the
  ideal point), dimension at
  a time.  If the legislators' party affiliations are availble in the
  \code{\link{rollcall}} object that was passed to \code{\link{ideal}},
  then legislators from the same party are plotted with a unique color.

  \code{type="votes"} produces a plot of classification rates for each
  roll call, by the percentage of legislators voting for the losing
  side.  The \code{x}-ordinate is jittered for clarity.
}
\value{
  After drawing plots on the current device, exits silently returning \code{invisible(NULL)}.
}

\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{predict.ideal}} \code{\link{ideal}} }

\examples{
data(s109)
f = system.file("extdata","id1.rda",package="pscl")
load(f)
phat <- predict(id1)
plot(phat,type="legis")
plot(phat,type="votes")
}
\keyword{hplot}% __ONLY ONE__ keyword per line
