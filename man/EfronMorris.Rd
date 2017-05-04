\name{EfronMorris}
\alias{EfronMorris}
\docType{data}
\title{Batting Averages for 18 major league baseball players, 1970}

\description{Batting averages for 18 major league baseball players,
  first 45 at bats of the 1970 season.}

\usage{data(EfronMorris)}

\format{
  \describe{
	\item{\code{name}}{character, name of player}
	\item{\code{team}}{character, team of player, abbreviated}
	\item{\code{league}}{character, National League or American League}
  	\item{\code{r}}{numeric, hits in 1st 45 at bats}
	\item{\code{y}}{numeric, \code{r}/45, batting average over 1st
  45 at bats}
  \item{\code{n}}{numeric, number of at bats, remainder of 1970
  season}
  \item{\code{p}}{numeric, batting average over remainder of 1970 season}
  }
}

\source{Efron, Bradley and Carl Morris. 1975. Data Analysis Using
Stein's Estimator and Its Generalizations.  \emph{Journal of the American
Statistical Association}. 70:311-319.
}
  
\examples{
data(EfronMorris)
attach(EfronMorris)
plot(p~y,
     xlim=range(p,y),
     ylim=range(p,y),
     xlab="Batting Average, 1st 45 at bats",
     ylab="Batting Average, Remainder of Season")
abline(0,1)
}
\keyword{datasets}
