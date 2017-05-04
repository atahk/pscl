\name{state.info}
\alias{state.info}
\docType{data}
\title{information about the American states needed for U.S. Congress}
\description{
  Numeric codes and names of 50 states and the District of Columbia,
  required to parse Keith Poole
  and Howard Rosenthal's collections of U.S. Congressional roll calls.
}
\usage{data(state.info)
}
\format{
  \describe{
    \item{\code{icpsr}}{integer, numeric code for state used by the
      Inter-university Consortium for Political and Social Research}
    \item{\code{state}}{character, name of state or \code{Washington
	D.C.}}
    \item{\code{year}}{numeric or \code{NA}, year of statehood}
    }
  }

\details{The function \code{\link{readKH}} converts the integer ICPSR codes
  into strings, via a table lookup in this data frame.  Another table
  lookup in \code{\link{state.abb}} provides the 2-letter abbreviation
  commonly used in identifying American legislators, e.g.,
  \code{KENNEDY, E (D-MA)}.} 
\seealso{\code{\link{state}}
}
\source{Various ICPSR codebooks. \url{http://www.icpsr.umich.edu}}
\keyword{datasets}
