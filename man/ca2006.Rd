\name{ca2006}
\alias{ca2006}
\docType{data}
\title{California Congressional Districts in 2006}
\description{
  Election returns and identifying information, California's 53
  congressional districts in the 2006 Congressional elections.
}
\usage{data(ca2006)}
\format{
  A data frame with 53 observations on the following 11 variables.
  \describe{
    \item{\code{district}}{numeric, number of Congressional district}
    \item{\code{D}}{numeric, number of votes for the Democratic candidate}
    \item{\code{R}}{numeric, votes for the Republican candidate}
    \item{\code{Other}}{numeric, votes for other candidates}
    \item{\code{IncParty}}{character, party of the incumbent (or retiring member),
      \code{D} or \code{R}}
    \item{\code{IncName}}{character, last name of the incumbent, character
      \code{NA} if no incumbent running}
    \item{\code{open}}{logical, \code{TRUE} if no incumbent running}
    \item{\code{contested}}{logical, \code{TRUE} if both major parties
      ran candidates}
    \item{\code{Bush2004}}{numeric, votes for George W. Bush (R) in the district in
      the 2004 presidential election}
    \item{\code{Kerry2004}}{numeric, votes for John Kerry (D) in 2004}
    \item{\code{Other2004}}{numeric votes for other candidates in 2004}
    \item{\code{Bush2000}}{numeric, votes for George W. Bush in 2000}
    \item{\code{Gore2000}}{numeric, votes for Al Gore (D) in 2000}
  }
}

\source{2006 data from the California Secretary of State's web site,
  \url{https://www.sos.ca.gov/elections/prior-elections/statewide-election-results/general-election-november-7-2006/statement-vote}. 
  %%Excel data at \url{http://elections.cdn.sos.ca.gov/sov/2006-general/congress.xls}. 
   2004 and 2000 presidential vote in congressional districts from the 2006 \emph{Almanac of American Politics}.

Thanks to Arthur Aguirre for the updated links, above.
}

\references{
  Michael Baraon and Richard E. Cohen.  2006. \emph{The Almanac of American
  Politics, 2006.} National Journal Group: Washington, D.C.
}
\examples{
data(ca2006)

## 2006 CA congressional vote against 2004 pvote
y <- ca2006$D/(ca2006$D+ca2006$R)
x <- ca2006$Kerry2004/(ca2006$Kerry2004+ca2006$Bush2004)

pch <- rep(19,length(y))
pch[ca2006$open] <- 1
col <- rep("black",length(y))
col[11] <- "red"    ## Pembo (R) loses to McNerney (D)
plot(y~x,pch=pch,
     col=col,
     xlim=range(x,y,na.rm=TRUE),
     ylim=range(x,y,na.rm=TRUE),
     xlab="Kerry Two-Party Vote, 2004",
     ylab="Democratic Two-Party Vote Share, 2006")
abline(0,1)
abline(h=.5,lty=2)
abline(v=.5,lty=2)
legend(x="topleft",
       bty="n",
       col=c("red","black","black"),
       pch=c(19,19,1),
       legend=c("Seat Changing Hands",
         "Seat Retained by Incumbent Party",
         "Open Seat (no incumbent running)")
       )
}
\keyword{datasets}
