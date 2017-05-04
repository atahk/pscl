\name{hitmiss}
\alias{hitmiss}
\alias{hitmiss.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Table of Actual Outcomes against Predicted Outcomes for discrete
  data models}

\description{ Cross-tabulations of actual outcomes against predicted
  outcomes for discrete data models, with summary statistics such as
  percent correctly predicted (PCP) under fitted and null models.  For
  models with binary responses (generalized linear models with
  \code{family=binomial}), the user can specific a classification
  threshold for the predicted probabilities.  }

\usage{
hitmiss(obj, digits = max(3, getOption("digits") - 3), ...)

\method{hitmiss}{glm}(obj,digits=max(3,getOption("digits")-3),
            ...,
            k=.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a fitted model object, such as a \code{glm} with
    \code{family=binomial}, a \code{polr} model for ordinal responses,
    or a \code{multinom} model for unordered/multinomial outcomes} 
  \item{digits}{number of digits to display in on-screen output}
  \item{\dots}{additional arguments passed to or from other functions}
  \item{k}{classification threshold for binary models}
}

\details{For models with binary responses, the user can specify a
  parameter 0 < \code{k} < 1; if the predicted probabilities exceed this
  threshold then the model is deemed to have predicted y=1, and
  otherwise to have predicted y=0.  Measures like percent correctly
  predicted are crude summaries of model fit; the cross-tabulation of
  actual against predicted is somewhat more informative, providing a
  little more insight as to where the model fits less well.
}

\value{
  For \code{hitmiss.glm}, a vector of length 3:
  \item{pcp}{Percent Correctly Predicted}
  \item{pcp0}{Percent Correctly Predicted among y=0}
  \item{pcp1}{Percent Correctly Predicted among y=1}
}

\author{Simon Jackman \email{jackman@stanford.edu}}

\note{To-do: The \code{glm} method should also handle binomial data presented
  as two-vector success/failures counts; and count data with
  \code{family=poisson}, the \code{glm.nb} models and \code{zeroinfl}
  and \code{hurdle} etc.   We should also make the output a class with
  prettier print methods, i.e., save the cross-tabulation in the
  returned object etc.
}

\seealso{\code{\link{pR2}} for pseudo r-squared;  \code{\link{predict}};
  \code{\link{extractAIC}}.  See also the \pkg{ROCR} package and the \code{lroc}  function in the \pkg{epicalc} package for ROC computations for assessing binary classifications.} 

\examples{
data(admit)
## ordered probit model
op1 <- MASS::polr(score ~ gre.quant + gre.verbal + ap + pt + female,
            Hess=TRUE,
            data=admit,
            method="probit")
hitmiss(op1)
}

\keyword{models}
