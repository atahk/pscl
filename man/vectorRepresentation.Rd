\name{vectorRepresentation}
\alias{vectorRepresentation}

\title{convert roll call matrix to series of vectors}

\description{
 Extract the information in a roll call matrix as a series of vectors
 with voting decision, a unique identifier for the legislator and a
 unique identifier for the roll call.
}

\usage{
vectorRepresentation(object, dropList = list(codes = c("missing", "notInLegis")))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{\link{rollcall}}}
  \item{dropList}{a \code{dropList}; see \code{\link{dropRollCall}}}
}

\details{ It is often the case that roll call matrices are sparse, say,
  when the roll call matrix has an \dQuote{overlapping generations}
  structure; e.g., consider forming data by pooling across a long
  temporal sequence of legislatures such that relatively few of the
  legislators in the data set actually vote on any given roll call.  In
  such a case, representing the data as a roll call matrix is not
  particularly helpful nor efficient, either for data summaries or modeling.
}

\value{A \code{\link{matrix}} with \code{z} rows, where \code{z} is the
  number of non-missing entries in \code{object$votes}, with
  \sQuote{missingness} defined by the \code{codes} component of the
  \code{dropList}.  The matrix has 3 columns:
  \item{vote}{the voting decision, either a \code{1} if the
    corresponding element of the roll call matrix \code{object$votes}
    is in the \code{yea} component of \code{object$codes}, or a \code{0}
    if the corresponding element of the roll call matrix is in the
    \code{nay} component of \code{object$codes}.  Non-missing entries of
    the roll call matrix are not stored.}
  \item{i}{the row of the roll call matrix \code{object$votes} that
    supplied the voting decision; i.e., a unique identifier for the
    legislator generating this \code{vote}}
  \item{j}{the column of the roll call matrix \code{object$votes} that
    supplied the \code{vote}; i.e., a unique identifier for the vote.}
}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\seealso{\code{\link{rollcall}}}
\examples{
data(s109)
y <- vectorRepresentation(s109)
apply(y,2,table,exclude=NULL)
}
\keyword{manip}
\keyword{utilities}
