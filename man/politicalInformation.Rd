\name{politicalInformation}
\alias{politicalInformation}
\docType{data}
\title{Interviewer ratings of respondent levels of political information}
\description{
  Interviewers administering the 2000 American National Election Studies
  assigned an ordinal rating to each respondent's "general
  level of information" about politics and public affairs.
}

\usage{data(politicalInformation)}

\format{
  A data frame with 1807 observations on the following 8 variables.
  \describe{
    \item{\code{y}}{interviewer rating, a factor with levels \code{Very Low} \code{Fairly Low} \code{Average} \code{Fairly High} \code{Very High}}
    \item{\code{collegeDegree}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{female}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{age}}{a numeric vector, respondent age in years}
    \item{\code{homeOwn}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{govt}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{length}}{a numeric vector, length of ANES pre-election
      interview in minutes}
    \item{\code{id}}{a factor, unique identifier for each interviewer}
  }
}

\details{
  Seven respondents have missing data on the ordinal interviewer rating.
   The covariates \code{age} and \code{length} also have some missing data.
}

\source{
  The National Election Studies (www.electionstudies.org). THE 2000
  NATIONAL ELECTION STUDY [dataset]. Ann Arbor, MI: University of
  Michigan, Center for Political Studies [producer and distributor].  
}

\references{
  Jackman, Simon. 2009.  \emph{Bayesian Analysis for the Social
  Sciences}.  Wiley: Hoboken, New Jersey.
}

\examples{
data(politicalInformation)

table(politicalInformation$y,exclude=NULL)

op <- MASS::polr(y ~ collegeDegree + female + log(age) + homeOwn + govt + log(length),
           data=politicalInformation,
           Hess=TRUE,
           method="probit")
}
\keyword{datasets}
