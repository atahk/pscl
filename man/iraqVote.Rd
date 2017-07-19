\name{iraqVote}
\alias{iraqVote}
\docType{data}
\title{
	U.S. Senate vote on the use of force against Iraq, 2002.	
}
\description{
	On October 11, 2002, the United States Senate
  voted 77-23 to authorize the use of military force against Iraq.  This data set lists the \dQuote{Ayes} and \dQuote{Nays} for each Senator and some covariates.
}
\usage{data(iraqVote)}
\format{
  A data frame with 100 observations on the following 6 variables.
  \describe{
    \item{\code{y}}{a numeric vector, the recorded vote (1 if Aye, 0 if Nay)}
    \item{\code{state.abb}}{two letter abbreviation for each state}
    \item{\code{name}}{senator name, party and state, e.g., \code{AKAKA (D HI)}}
    \item{\code{rep}}{logical, \code{TRUE} for Republican senators}
    \item{\code{state.name}}{name of state}
    \item{\code{gorevote}}{numeric, the vote share recorded by Al Gore in the corresponding state in the 2000 Presidential election}
  }
}
\details{The only Republican to vote against the resolution was Lincoln
  Chafee (Rhode Island); Democrats split 29-22 in favor of the
  resolution.
}
\source{
	Keith Poole, 107th Senate Roll Call Data. \url{https://voteview.com/static/data/out/votes/S107_votes.ord}  The Iraq vote is vote number 617.

	David Leip's Atlas of U.S. Presidential Elections.  \url{http://uselectionatlas.org} 
}
\references{
	Jackman, Simon. 2009.  \emph{Bayesian Analysis for the Social Sciences}.  Wiley: Chichester.  Example 8.3.
}
\examples{
data(iraqVote)
## probit model
glm1 <- glm(y ~ gorevote + rep,
            data=iraqVote,
            family=binomial(link=probit))
}
\keyword{datasets}
