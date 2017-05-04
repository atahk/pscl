\name{sc9497}
\docType{data}
\alias{sc9497}
\title{votes from the United States Supreme Court, from 1994-1997}
\description{
  This data set provides information on the United States Supreme Court
  from 1994-1997.  Votes included are non-unanimous.
}

\usage{data(sc9497)}
\format{A list containing the elements:
  \describe{
    \item{votes}{a matrix of the votes, 0=Nay, 1=Yea,
    \code{NA}=Abstained or missing data.  The matrix columns are labeled with
    \code{vote.names} and the rows are labeled with \code{legis.names}.}
    \item{legis.names}{a vector of the names of the nine Justices
    sitting on the court at this time.}
    \item{party}{NULL; exists for consistency with House and Senate data sets.}
    \item{state}{NULL; exists for consistency with House and Senate data
    sets.}
    \item{district}{NULL; exists for consistency with House data sets.}
    \item{id}{NULL; exists for consistency with House and Senate data sets.}
    \item{vote.names}{a vector of strings numbering the cases simply to
      distinguish them from one another.}
    \item{desc}{a description of the data set.}
  }
}

\source{Harold J. Spaeth (1999). \emph{United States Supreme Court Judicial
Database, 1953-1997 Terms}.  Ninth edition.  Inter-university Consortium for 
                        Political and Social Research.  Ann Arbor,
			Michigan.  \url{http://www.icpsr.umich.edu/}}

\keyword{datasets}
