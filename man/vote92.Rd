\name{vote92}
\alias{vote92}
\docType{data}
\title{
	Reports of voting in the 1992 U.S. Presidential election.
}
\description{
	Survey data containing self-reports of vote choice in the 1992 U.S. Presidential election, with numerous covariates, from the 1992 American National Election Studies.
}
\usage{data(vote92)}
\format{
  A data frame with 909 observations on the following 10 variables.
  \describe{
    \item{\code{vote}}{a factor with levels \code{Perot} \code{Clinton} \code{Bush}}
    \item{\code{dem}}{a numeric vector, 1 if the respondent reports identifying with the Democratic party, 0 otherwise.}
    \item{\code{rep}}{a numeric vector, 1 if the respondent reports identifying with the Republican party, 0 otherwise}
    \item{\code{female}}{a numeric vector, 1 if the respondent is female, 0 otherwise}
    \item{\code{persfinance}}{a numeric vector, -1 if the respondent reports that their personal financial situation has gotten worse over the last 12 months, 0 for no change, 1 if better}
    \item{\code{natlecon}}{a numeric vector, -1 if the respondent reports that national economic conditions have gotten worse over the last 12 months, 0 for no change, 1 if better}
    \item{\code{clintondis}}{a numeric vector, squared difference between respondent's self-placement on a scale measure of political ideology and the respondent's placement of the Democratic candidate, Bill Clinton}
    \item{\code{bushdis}}{a numeric vector, squared ideological distance of the respondent from the Republican candidate, President George H.W. Bush}
    \item{\code{perotdis}}{a numeric vector, squared ideological distance of the respondent from the Reform Party candidate, Ross Perot}
  }
}
\details{These data are unweighted.  Refer to the original data source for weights that purport to correct for non-representativeness and non-response.}
\source{
	Alvarez, R. Michael and Jonathan Nagler. 1995. Economics, issues and the Perot candidacy: Voter choice in the 1992 Presidential election.  \emph{American Journal of Political Science}. 39:714-44.

	Miller, Warren E., Donald R. Kinder, Steven J. Rosenstone and the National Election Studies.  1999.  \emph{National Election Studies, 1992: Pre-/Post-Election Study}.  Center for Political Studies, University of Michigan: Ann Arbor, Michigan.

	 Inter-University Consortium for Political and Social Research.  Study Number 1112.  \url{http://dx.doi.org/10.3886/ICPSR01112}.   	
}
\references{
	Jackman, Simon. 2009.  \emph{Bayesian Analysis for the Social Sciences}.
	Wiley: Hoboken, New Jersey.  Examples 8.7 and 8.8.
}
\examples{
data(vote92)
summary(vote92)
}
\keyword{datasets}
