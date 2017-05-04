\name{hurdle}
\alias{hurdle}
\alias{print.hurdle}

\title{Hurdle Models for Count Data Regression}

\description{
  Fit hurdle regression models for count data via maximum likelihood. 
}

\usage{
hurdle(formula, data, subset, na.action, weights, offset,
  dist = c("poisson", "negbin", "geometric"),
  zero.dist = c("binomial", "poisson", "negbin", "geometric"),
  link = c("logit", "probit", "cloglog", "cauchit", "log"),
  control = hurdle.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model, see details.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor of the count model. See below for more
    information on offsets.}
  \item{dist}{character specification of count model family.}
  \item{zero.dist}{character specification of the zero hurdle model family.}
  \item{link}{character specification of link function in the binomial
    zero hurdle (only used if \code{zero.dist = "binomial"}.}
  \item{control}{a list of control arguments specified via
    \code{\link{hurdle.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.}
  \item{\dots}{arguments passed to \code{\link{hurdle.control}} in the
    default setup.}
}

\details{
  Hurdle count models are two-component models with a truncated count
  component for positive counts and a hurdle component that models the
  zero counts. Thus, unlike zero-inflation models, there are \emph{not} two
  sources of zeros: the count model is only employed if the hurdle for
  modeling the occurence of zeros is exceeded. The count model is typically
  a truncated Poisson or negative binomial regression (with log link).
  The geometric distribution is a special case of the negative binomial with
  size parameter equal to 1. For modeling the hurdle (occurence of positive
  counts) either a binomial model can be employed or a censored count distribution.
  Binomial logit and censored geometric models as the hurdle part both lead to 
  the same likelihood function and thus to the same coefficient estimates.
  A censored negative binomial model for the zero hurdle is only identified
  if there is at least one non-constant regressor with (true) coefficient different
  from zero (and if all coefficients are close to zero the model can be poorly
  conditioned).
  
  The \code{formula} can be used to specify both components of the model:
  If a \code{formula} of type \code{y ~ x1 + x2} is supplied, then the same
  regressors are employed in both components. This is equivalent to
  \code{y ~ x1 + x2 | x1 + x2}. Of course, a different set of regressors
  could be specified for the zero hurdle component, e.g.,
  \code{y ~ x1 + x2 | z1 + z2 + z3} giving the count data model \code{y ~ x1 + x2}
  conditional on (\code{|}) the zero hurdle model \code{y ~ z1 + z2 + z3}.
  
  Offsets can be specified in both parts of the model pertaining to count and
  zero hurdle model: \code{y ~ x1 + offset(x2) | z1 + z2 + offset(z3)}, where
  \code{x2} is used as an offset (i.e., with coefficient fixed to 1) in the
  count part and \code{z3} analogously in the zero hurdle part. By the rule
  stated above \code{y ~ x1 + offset(x2)} is expanded to
  \code{y ~ x1 + offset(x2) | x1 + offset(x2)}. Instead of using the
  \code{offset()} wrapper within the \code{formula}, the \code{offset} argument
  can also be employed which sets an offset only for the count model. Thus,
  \code{formula = y ~ x1} and \code{offset = x2} is equivalent to
  \code{formula = y ~ x1 + offset(x2) | x1}.
  
  All parameters are estimated by maximum likelihood using \code{\link[stats]{optim}},
  with control options set in \code{\link{hurdle.control}}.
  Starting values can be supplied, otherwise they are estimated by \code{\link[stats]{glm.fit}}
  (the default). By default, the two components of the model are estimated separately
  using two \code{optim} calls. Standard errors are derived numerically using
  the Hessian matrix returned by \code{\link[stats]{optim}}. See
  \code{\link{hurdle.control}} for details.
  
  The returned fitted model object is of class \code{"hurdle"} and is similar
  to fitted \code{"glm"} objects. For elements such as \code{"coefficients"} or
  \code{"terms"} a list is returned with elements for the zero and count components,
  respectively. For details see below.
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"hurdle"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{fitted}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.matrix}}. See \code{\link{predict.hurdle}} for more details
  on all methods.
}

\value{
  An object of class \code{"hurdle"}, i.e., a list with components including
  \item{coefficients}{a list with elements \code{"count"} and \code{"zero"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{a list (of lists) with the output(s) from the \code{optim} call(s) for
    minimizing the negative log-likelihood(s),}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call(s),}
  \item{weights}{the case weights used,}
  \item{offset}{a list with elements \code{"count"} and \code{"zero"}
    containing the offset vectors (if any) from the respective models,}
  \item{n}{number of observations (with weights > 0),}
  \item{df.null}{residual degrees of freedom for the null model (= \code{n - 2}),}
  \item{df.residual}{residual degrees of freedom for fitted model,}
  \item{terms}{a list with elements \code{"count"}, \code{"zero"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{theta}{estimate of the additional \eqn{\theta}{theta} parameter of the
    negative binomial model(s) (if negative binomial component is used),}
  \item{SE.logtheta}{standard error(s) for \eqn{\log(\theta)}{log(theta)},}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of all coefficients in the model (derived from the
    Hessian of the \code{optim} output(s)),}
  \item{dist}{a list with elements \code{"count"} and \code{"zero"} with character
    strings describing the respective distributions used,}
  \item{link}{character string describing the link if a binomial zero hurdle model
    is used,}
  \item{linkinv}{the inverse link function corresponding to \code{link},}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"count"} and \code{"zero"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response count vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"count"} and \code{"zero"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\references{
Cameron, A. Colin and Pravin K. Trivedi. 1998. \emph{Regression Analysis of Count 
Data}. New York: Cambridge University Press.

Cameron, A. Colin and Pravin K. Trivedi 2005. \emph{Microeconometrics: Methods and Applications}.
Cambridge: Cambridge University Press.

Mullahy, J. 1986. Specification and Testing of Some Modified Count Data Models.
  \emph{Journal of Econometrics}. \bold{33}:341--365.

Zeileis, Achim, Christian Kleiber and Simon Jackman 2008.
\dQuote{Regression Models for Count Data in R.} 
\emph{Journal of Statistical Software}, \bold{27}(8).
URL \url{http://www.jstatsoft.org/v27/i08/}.
}

\author{Achim Zeileis <Achim.Zeileis@R-project.org>}

\seealso{\code{\link{hurdle.control}}, \code{\link[stats]{glm}},
  \code{\link[stats]{glm.fit}}, \code{\link[MASS]{glm.nb}},
  \code{\link{zeroinfl}}
}

\examples{
## data
data("bioChemists", package = "pscl")

## logit-poisson
## "art ~ ." is the same as "art ~ . | .", i.e.
## "art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment"
fm_hp1 <- hurdle(art ~ ., data = bioChemists)
summary(fm_hp1)

## geometric-poisson
fm_hp2 <- hurdle(art ~ ., data = bioChemists, zero = "geometric")
summary(fm_hp2)

## logit and geometric model are equivalent
coef(fm_hp1, model = "zero") - coef(fm_hp2, model = "zero")

## logit-negbin
fm_hnb1 <- hurdle(art ~ ., data = bioChemists, dist = "negbin")
summary(fm_hnb1)

## negbin-negbin
## (poorly conditioned zero hurdle, note the standard errors)
fm_hnb2 <- hurdle(art ~ ., data = bioChemists, dist = "negbin", zero = "negbin")
summary(fm_hnb2)
}

\keyword{regression}
