\name{seatsVotes}
\alias{seatsVotes}
%\alias{print.seatsVotes}
%\alias{summary.seatsVotes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A class for creating seats-votes curves}
\description{
  Convert a vector of vote shares into a seats-vote curve object,
  providing estimates of partisan bias.
}
\usage{
seatsVotes(x, desc = NULL, method = "uniformSwing")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of vote shares for a specific party (either
    proportions or percentages)}
  \item{desc}{descriptive text}
  \item{method}{how to simulate a seats-vote curve; the only supported
    method at this stage is \code{uniformSwing}.}
}

\details{Simulation methods are required to induce a seats-votes curve
  given a vector of vote shares from one election.  The uniform swing
  method simply slides the empirical distribution function of the vote
  shares \dQuote{up} and \dQuote{down}, computing the proportion of
  the vote shares that lie above .5 (by construction, the winning
  percentage in a two-party election) for each new location of the
  vector of vote shares.  That is, as the empirical CDF of the
  observed vote shares slides up or down, more or less seats cross the
  .5 threshold.  A seats-votes curve is formed by plotting the seat
  share above .5 as a function of the average district-level vote
  share (a weakly monotone function, since the empirical CDF
  constitutes a set of sufficient statistics for this problem).  The
  simulation is run so as to ensure that average district-level vote
  shares range between 0 and 1.

  The extent to which the seats-votes curve departs from symmetry is
  known as bias.  More specifically, the vertical displacement of the
  seats-votes curve from .5 when average district-level vote share is .5
  is conventionally reported as an estimate of the bias of the electoral
  system.  

  Different methods produce different estimates of seats-votes curves
  and summary estimands such as bias.  The
  uniform swing method is completely deterministic and does not
  produce any uncertainty assessment (e.g., confidence intervals etc).
}

\value{
  An object of class \code{seatsVotes}, with components
  \item{s}{Estimated seat shares over the range of simulated average,
    district-level vote shares}
  \item{v}{Simulated average district-level vote shares}
  \item{x}{observed seat shares, with missing data removed}
  \item{desc}{user-supplied descriptive character string}
  \item{call}{a list of class \code{\link{call}}, the call to the function}
}

\references{
   Tufte, Edward R. 1973. The Relationship Between Seats and Votes in
   Two-Party Systems. \emph{American Political Science
     Review}. 67(2):540-554.

   Gelman, Andrew and Gary King.  1990.  Estimating the Consequences of
   Electoral Redistrictring.  \emph{Journal of the American Statistical
     Association}. 85:274-282.

   Jackman, Simon. 1994.  Measuring Electoral Bias: Australia,
   1949-93. \emph{British Journal of Political Science}. 24(3):319-357.
}

\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\note{Additional methods to come later.}

\seealso{\code{\link{plot.seatsVotes}} for plotting methods.}

\examples{
data(ca2006)
x <- ca2006$D/(ca2006$D+ca2006$R)
sv <- seatsVotes(x,
                 desc="Democratic Vote Shares, California 2006 congressional elections")
}

\keyword{misc}
