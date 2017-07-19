\name{predprob}
\alias{predprob}

\title{compute predicted probabilities from fitted models}

\description{
  Compute predicted probabilities from fitted models, optionally at new
  covariate values.
}

\usage{
predprob(obj, \dots)
}

\arguments{
  \item{obj}{fitted model object}
  \item{\dots}{other arguments}
}

\details{
  See documentation for specific methods.
}

\value{
  A matrix of predicted probabilities, each row a vector of
  predicted probabilities over the range of responses seen in the data
  (i.e., \code{min(y):max(y)}), conditional on the values of covariates. 
}

\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}

\seealso{\code{\link{predprob.glm}}, \code{\link{predprob.zeroinfl}}}

\examples{
data("bioChemists")
zip <- zeroinfl(art ~ . | ., data = bioChemists, EM = TRUE)
phat <- predprob(zip)

newdata <- expand.grid(list(fem="Men",mar="Married",
                            kid5=1,phd=3.103,
                            ment=0:77))
phat <- predprob(zip, newdata = newdata)
}
\keyword{regression}
\keyword{models}
