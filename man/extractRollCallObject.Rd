\name{extractRollCallObject}
\alias{extractRollCallObject}
\title{return the roll call object used in fitting an ideal model}
\description{
  Given a fitted model of class \code{\link{ideal}}, return the
  \code{\link{rollcall}} object that was used in the model fitting
  (i.e., apply all subseting and recoding implied by the \code{dropList}
  passed to \code{\link{ideal}}).
}
\usage{
extractRollCallObject(object)
}
\arguments{
  \item{object}{an object of class \code{\link{ideal}}}
}
\details{
  This function is used by many post-estimation commands that operate on
  objects of class \code{\link{ideal}}.  The function inspects the
  \code{call} attribute of the 
  \code{\link{ideal}} object, extracting the name of the
  \code{\link{rollcall}} object and the \code{dropList}, then hands them
  over to \code{\link{dropRollCall}}.
}
\value{
  An object of class \code{\link{rollcall}}
}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\seealso{\code{\link{rollcall}}; see \code{\link{dropRollCall}} for
  details on the form of a \code{dropList}.}
\examples{
data(s109)
f = system.file("extdata","id1.rda",package="pscl")
load(f)
tmp <- extractRollCallObject(id1)
summary(tmp)
v <- convertCodes(tmp)         ## roll call matrix per se  
}
\keyword{models}
