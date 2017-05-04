\name{idealToMCMC}
\alias{idealToMCMC}
\title{convert an object of class ideal to a coda MCMC object}
\description{
  Converts the x element of an \code{ideal} object to an MCMC object,
  as used in the \pkg{coda} package.
}

\usage{
idealToMCMC(object, burnin=NULL)
}

\arguments{
  \item{object}{an object of class \code{\link{ideal}}.}

 \item{burnin}{of the recorded MCMC samples, how many to discard as
    burnin?  Default is \code{NULL}, in which case the value of
    \code{burnin} in the \code{\link{ideal}} object is used.}
}

\value{A \code{\link[coda:mcmc]{mcmc}} object as used by the
  \pkg{coda} package, starting at
  iteration \code{start}, drawn from the \code{x} component of the
  \code{\link{ideal}} object.}

\note{When specifying a value of \code{burnin} different from that used
  in fitting the \code{\link{ideal}} object, note a distinction
  between the iteration numbers of the stored iterations, and the
  number of stored iterations.  That is, the \code{n}-th iteration
  stored in an \code{\link{ideal}} object will not be iteration
  \code{n} if the user specified \code{thin>1} in the call to
  \code{\link{ideal}}.  Here, iterations are tagged with their
  iteration number.  Thus, if the user called \code{\link{ideal}} with
  \code{thin=10} and \code{burnin=100} then the stored iterations are
  numbered \code{100, 110, 120, ...}.  Any future subsetting via a
  \code{burnin} refers to this iteration number.}

\seealso{\code{\link{ideal}}, \code{\link[coda:mcmc]{mcmc}}}

\examples{
data(s109)
f = system.file("extdata",package="pscl","id1.rda")
load(f)
id1coda <- idealToMCMC(id1)
summary(id1coda)  
}
\keyword{classes}

