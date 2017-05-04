\name{s109}
\alias{s109}
\docType{data}
\title{rollcall object, 109th U.S. Senate (2005-06).}

\description{A sample rollcall object, generated using a
  collection of the rollcalls of the 109th U.S. Senate (2005-2006).}

\usage{data(s109)}

\format{
  A \code{\link{rollcall}} object containing the recorded votes of the 109th
  U.S. Senate, plus information identifying the legislators and the
  rollcalls.
}

\details{Note the coding scheme used by Poole and Rosenthal; Yea (1),
  Nay (6) etc.}

\source{
  Keith Poole's web site:
  \url{http://voteview.com/senate109.htm}

  Originally scraped from the Senate's web site by Jeff Lewis (UCLA).
  
  Information identifying the votes is available at
  \url{ftp://voteview.com/dtaord/s109desc_final.csv}
}

\examples{
require(pscl)
data(s109)
is(s109,"rollcall")    ## TRUE
s109                   ## print method for class rollcall
summary(s109)          ## summary method
summary(s109,verbose=TRUE)
\dontrun{
## how s109 was created
require(pscl)
s109 <- readKH("ftp://voteview.com/dtaord/sen109kh.ord",
               desc="109th U.S. Senate",
               debug=TRUE)
url <- "ftp://voteview.com/dtaord/s109desc_final.csv"

s109$vote.data <- data.frame(read.csv(file=url,header=TRUE))
s109$vote.data$date <- as.Date(s109$vote.data$date,
                               format="%Y-%m-%d")
dimnames(s109$votes)[[2]] <- paste(s109$vote.data$session,
                                   s109$vote.data$number,sep="-")
}
}
\keyword{datasets}
