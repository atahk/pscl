\name{vuong}
\alias{vuong}
\title{Vuong's non-nested hypothesis test}
\description{
  Compares two models fit to the same data that do not nest via Vuong's
  non-nested test.
}

\usage{
vuong(m1, m2, digits = getOption("digits"))
}

\arguments{
  \item{m1}{model 1, an object inheriting from class \code{glm},
    \code{negbin} or \code{zeroinfl}}
  \item{m2}{model 2, as for model 1}
  \item{digits}{significant digits in printed result}
}
\details{
  The Vuong non-nested test is based on a comparison of the predicted
  probabilities of two models that do not nest.  Examples include
  comparisons of zero-inflated count models with their non-zero-inflated
  analogs (e.g., zero-inflated Poisson versus ordinary Poisson, or
  zero-inflated negative-binomial versus ordinary negative-binomial).  A
  large, positive test statistic provides evidence of
  the superiority of model 1 over model 2, while a large, negative test statistic is
  evidence of the superiority of model 2 over model 1.  Under the null
  that the models are indistinguishable, the test
  statistic is asymptotically distributed standard normal.
  
  Let \eqn{p_i = \hat{Pr}(y_i | M_1)} be the predicted probabilities from model 1, evaluated conditional on the estimated MLEs.   Let \eqn{q_i} be the corresponding probabilities from model 2.   Then the Vuong statistic is \eqn{\sqrt{N} \bar{m}/s_m} where \eqn{m_i = log(p_i) - log(q_i)} and \eqn{s_m} is the sample standard deviation of \eqn{m_i}.
  
  Two finite sample corrections are often considered, based on the Akaike (AIC) and Schwarz (BIC) penalty terms, based on the complexity of the two models. These corrections sometimes generate conflicting conclusions. 

  The function will fail if the models do not contain identical values
  in their respective components named \code{y} (the value of the
  response being modeled).
  
  
}
\value{
  nothing returned, prints 3 test statistics and \eqn{p} values and
  exits silently.
}
\references{Vuong, Q.H. 1989. Likelihood ratio tests for model
  selection and non-nested hypotheses. \emph{Econometrica}. 57:307-333.}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}

\examples{
\dontrun{
data("bioChemists")
## compare Poisson GLM and ZIP
glm1 <- glm(art ~ ., data = bioChemists, family = poisson)
zip <- zeroinfl(art ~ . | ., data = bioChemists, EM = TRUE)
vuong(glm1, zip)


## compare negbin with zero-inflated negbin
nb1 <- MASS::glm.nb(art ~ ., data=bioChemists)
zinb <- zeroinfl(art ~ . | ., data = bioChemists, dist = "negbin", EM = TRUE)
vuong(nb1, zinb)
}
}
\keyword{models}
