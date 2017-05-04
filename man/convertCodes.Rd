\name{convertCodes}
\alias{convertCodes}
\title{convert entries in a rollcall matrix to binary form}

\description{
  Convert roll call matrix to binary form using encoding information.
}

\usage{
convertCodes(object, codes = object$codes)
}

\arguments{
  \item{object}{\code{\link{rollcall}} object}
  \item{codes}{list, mapping entries in the \code{votes} component of
    \code{rollcall} object to 0 (\sQuote{Nay}), 1 (\sQuote{Yea}) and
    \code{NA} (missing, abstentions, etc). Defaults to the \code{codes}
    component of the \code{rollcall} object.} 
}

\details{
  See \code{\link{rollcall}} for details on the form of the
  \code{codes} list.
}
\value{a \code{\link{matrix}} with dimensions equal to the dimensions of
  the \code{votes} component of the \code{rollcall} object.} 
\author{Simon Jackman \email{jackman@stanford.edu}}
\note{Any entries in the \code{votes} matrix that can
  not be mapped into \code{c(0,1,NA)} using the information in
  \code{codes} are mapped to \code{NA}, with an informative message
  sent to the console.
}
\seealso{\code{\link{rollcall}}}
\examples{
  data(s109)
  mat <- convertCodes(s109)
  table(mat,exclude=NULL)
}  
\keyword{manip}
