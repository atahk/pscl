\name{nj07}
\alias{nj07}
\docType{data}
\title{rollcall object, National Journal key votes of 2007}

\description{A rollcall object containing 99 rollcalls from the 2nd session of the 110th U.S. Senate,  
  designated by \emph{National Journal} as the "key votes" of 2007.  These data were used to by \emph{National Journal} to rate (then Senator) Barack Obama was the "most liberal senator" in 2007.} 

\usage{data(nj07)}

\format{
  A \code{\link{rollcall}} object containing the recorded votes,
  plus information identifying the legislators and the
  rollcalls.
}

\details{Note the coding scheme used by Poole and Rosenthal; Yea (1,2,3),
  Nay (4,5,6) etc.}

\source{
  Keith Poole's web site:
  \url{https://legacy.voteview.com/senate110.htm}

  Originally scraped from the Senate's web site by Jeff Lewis.

  Josh Clinton compiled the list of \emph{National Journal} key votes.  
}

\references{
	Clinton, Joshua and Simon Jackman.  2009. To Simulate or NOMINATE?  
	\emph{Legislative Studies Quarterly}. V34(4):593-621.

	Jackman, Simon.  2009. \emph{Bayesian Analysis for the Social Sciences}. Wiley: Hoboken, New Jersey. Example 9.2. 
}

\examples{
require(pscl)
data(nj07)
is(nj07,"rollcall")    ## TRUE
nj07                   ## print method for class rollcall
names(nj07)
names(nj07$vote.data)
table(nj07$vote.data$policyArea)
summary(nj07)          ## summary method
summary(nj07,verbose=TRUE)
}
\keyword{datasets}
