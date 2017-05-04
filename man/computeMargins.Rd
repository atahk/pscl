\name{computeMargins}
\alias{computeMargins}
\title{add information about voting outcomes to a rollcall
  object}

\description{
  Add summaries of each roll call vote to a \code{\link{rollcall}} object.
}

\usage{
computeMargins(object, dropList = NULL)
}

\arguments{
  \item{object}{an object of class \code{\link{rollcall}}}
  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}})
      listing voting decisions, legislators and/or votes to be dropped
      from the analysis; see \code{\link{dropRollCall}} for details.}
}

\details{ The subseting implied by the \code{dropList} is first applied
  to the \code{\link{rollcall}} object, via \code{\link{dropRollCall}}.
  Then, for each remaining roll call vote, the number of legislators
  voting \dQuote{Yea}, \dQuote{Nay}, and not voting are computed, using
  the encoding information in the \code{codes} component of the
  \code{\link{rollcall}} object via the \code{\link{convertCodes}}
  function.  The matrix of vote counts are added to the \code{\link{rollcall}}
  object as a component \code{voteMargins}.
}

\value{
  An object of class \code{\link{rollcall}}, with a component
  \code{voteMargins} that is a matrix with four columns: 
  \item{Yea}{number of legislators voting \dQuote{Yea}}
  \item{Nay}{number of legislators voting \dQuote{Nay}}
  \item{NA}{number of legislators not voting \dQuote{Nay}}
  \item{Min}{the number of legislators voting on the losing side of the
  roll call}
}

\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{dropRollCall}} on specifying a \code{dropList}.
  The vote-specific marginals produced by theis function are used by as
  \code{\link{dropRollCall}}, \code{\link{summary.ideal}} and
  \code{\link{predict.ideal}}.}

\examples{
data(s109)
tmp <- computeMargins(s109)
dim(tmp$voteMargins)   ## 645 by 4

tmp <- computeMargins(s109,
                     dropList=list(codes="notInLegis",lop=0))
dim(tmp$voteMargins)   ## 544 by 4
}

\keyword{manip} 
