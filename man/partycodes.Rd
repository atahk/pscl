\name{partycodes}
\alias{partycodes}
\docType{data}
\title{political parties appearing in the U.S. Congress}

\description{
  Numeric codes and names of 85 political parties appearing in Poole
  and Rosenthal's collection of U.S. Congressional roll calls.
}

\usage{data(partycodes)}

\format{
  \describe{
    \item{\code{code}}{integer, numeric code for legislator appearing in Poole and
      Rosenthal rollcall data files}
    \item{\code{party}}{character, name of party}
    }
  }

\details{The function \code{\link{readKH}} converts the integer codes
  into strings, via a table lookup in this data frame.}
\seealso{\code{\link{readKH}}
}
\source{Keith Poole's website: \url{https://legacy.voteview.com/PARTY3.HTM}}
\keyword{datasets}
