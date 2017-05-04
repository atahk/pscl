\name{admit}
\alias{admit}
\docType{data}
\title{Applications to a Political Science PhD Program}
\description{
  Ordinal ratings (faculty evaluations) of applicants to a Political Science PhD
  Program.}
\usage{data(admit)}
\format{
  A data frame with 106 observations on the following 6 variables.
  \describe{
    \item{\code{score}}{an ordered factor with levels \code{1} < \code{2} < \code{3} < \code{4} < \code{5}}
    \item{\code{gre.quant}}{applicant's score on the quantitative
      section of the GRE; the maximum score is 800}
    \item{\code{gre.verbal}}{applicant's score on the verbal section of
      the GRE; the maximum score is 800}
    \item{\code{ap}}{1 if the applicant indicated an interest in
      American politics; 0 otherwise}
    \item{\code{pt}}{1 if the applicant indicated an interest in
      Political Theory; 0 otherwise}
    \item{\code{female}}{1 for female applicants; 0 otherwise}
  }
}

\references{Jackman, Simon.  2004.  "What Do We Learn From Graduate
  Admissions Committees?: A Multiple-Rater, Latent Variable Model, with
  Incomplete Discrete and Continuous Indicators."  \emph{Political
  Analysis}.  12(4):400-424.  }

\examples{
data(admit)
summary(admit)
## ordered probit model
op1 <- MASS::polr(score ~ gre.quant + gre.verbal + ap + pt + female,
            Hess=TRUE,
            data=admit,
            method="probit")
summary(op1)
hitmiss(op1)
logLik(op1)
pR2(op1)   
}
\keyword{datasets}
