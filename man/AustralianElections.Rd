\name{AustralianElections}
\alias{AustralianElections}
\docType{data}
\title{elections to Australian House of
  Representatives, 1949-2016}

\description{
  Aggregate data on the 24 elections to Australia's House of
  Representatives, 1949 to 2016. 
}

\usage{data(AustralianElections)}

\format{
 A data frame with the following variables:	
  \describe{
    \item{\code{date}}{date of election, stored using the
      \code{\link{Date}} class}
    \item{\code{Seats}}{numeric, number of seats in the House of
      Representatives} 
    \item{\code{Uncontested}}{numeric, number of uncontested seats}
    \item{\code{ALPSeats}}{numeric, number of seats won by the
      Australian Labor Party}
    \item{\code{LPSeats}}{numeric, number of seats won by the Liberal Party} 
    \item{\code{NPSeats}}{numeric, number of seats won by the National
      Party (previously known as the Country Party)}
    \item{\code{OtherSeats}}{numeric, number of seats won by other
      parties and/or independent candidates}
    \item{\code{ALP}}{numeric, percentage of first preference votes cast for
      Australian Labor Party candidates}
    \item{\code{ALP2PP}}{numeric, percentage of the two-party preferred
      vote won by Australian Labor Party candidates}
    \item{\code{LP}}{numeric, percent of first preference votes cast for
      Liberal Party candidates}
    \item{\code{NP}}{numeric, percent of first preference votes cast for
      National Party (Country Party) candidates}
    \item{\code{DLP}}{numeric, percent of first preference votes cast
      for Democratic Labor Party candidates}
    \item{\code{Dem}}{numeric, percent of first preference votes cast
      for Australian Democrat candidates}
    \item{\code{Green}}{numeric, percent of first preference votes cast
      for Green Party candidates}
    \item{\code{Hanson}}{numeric, percent of first preference votes cast
      for candidates from Pauline Hanson's One Nation party}
    \item{\code{Com}}{numeric, percent of first preference votes cast
      for Communist Party candidates}
    \item{\code{AP}}{numeric, percent of first preference votes cast for
      Australia Party candidates}
    \item{\code{Informal}}{numeric, percent of ballots cast that are
      spoiled, blank, or otherwise uncountable (usually because of
      errors in enumerating preferences)}
    \item{\code{Turnout}}{numeric, percent of enrolled voters recorded
      as having turned out to vote (Australia has compulsory voting)}
  }
}

\note{The Liberal National Party of Queensland formed in 2008 after a merger of the Liberal Party and the National Party. In all elections following 2008, they have been categorised under \code{LP}.}

\source{Australian Electoral Commission. \url{http://www.aec.gov.au}.
}

\references{Jackman, Simon. 2009. \emph{Bayesian Analysis for the Social Sciences}.  Wiley: Hoboken, New Jersey.  Example 3.5.
} 
 
\examples{
data(AustralianElections)
attach(AustralianElections)
alpSeatShare <- ALPSeats/Seats
alpVoteShare <- ALP2PP/100

## log-odds transforms
x <- log(alpVoteShare/(1-alpVoteShare))
y <- log(alpSeatShare/(1-alpSeatShare))

ols <- lm(y~x)   ## Tufte-style seats-votes regression

xseq <- seq(-4.5,4.5,length=500)
yhat <- coef(ols)[1] + coef(ols)[2]*xseq
yhat <- exp(yhat)/(1+exp(yhat))
xseq <- exp(xseq)/(1+exp(xseq))

## seats vote curve
plot(x=alpVoteShare,
     y=alpSeatShare,
     xlab="ALP Vote Share",
     ylab="ALP Seat Share")
lines(xseq,yhat,lwd=2)
abline(h=.5,lty=2)
abline(v=.5,lty=2)
}
\keyword{datasets}
