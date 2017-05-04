\name{betaHPD}
\alias{betaHPD}

\title{compute and optionally plot beta HDRs}

\description{Compute and optionally plot highest density regions for the Beta
  distribution.}

\usage{
   betaHPD(alpha,beta,p=.95,plot=FALSE,xlim=NULL,debug=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{scalar, first shape parameter of the Beta density.  Must be
    greater than 1, see details}
  \item{beta}{scalar, second shape parameter of the Beta density.  Must be
    greater than 1, see details}
  \item{p}{scalar, content of HPD, must lie between 0 and 1}
  \item{plot}{logical flag, if \code{TRUE} then plot the density and
    show the HDR}
  \item{xlim}{numeric vector of length 2, the limits of the density's
   support to show when plotting; the default is \code{NULL}, in which
   case the function will confine plotting to where the density is
   non-neglible} 
  \item{debug}{logical flag, if \code{TRUE} produce messages to the
   console}
}
\details{The Beta density arises frequently in Bayesian models of
  binary events, rates, and proportions, which take on values in the
  open unit interval.  For instance, the Beta density is a conjugate prior
  for the unknown success probability in binomial trials.  With shape
  parameters \eqn{\alpha > 1} and \eqn{\beta > 1}, the Beta density is
  unimodal.

  In general, suppose \eqn{\theta \in \Theta \subseteq R^k}
  is a random variable with density \eqn{f(\theta)}.  A highest
  density region (HDR) of \eqn{f(\theta)} with content \eqn{p \in
  (0,1]} is a set \eqn{\mathcal{Q} \subseteq \Theta} with the
  following properties:
  \deqn{\int_\mathcal{Q} f(\theta) d\theta = p}
  and
  \deqn{f(\theta) > f(\theta^*) \, \forall\
    \theta \in \mathcal{Q},
    \theta^* \not\in \mathcal{Q}.}
  For a unimodal
  Beta density (the class of Beta densities handled by this function),
  a HDR of content \eqn{0 < p < 1} is simply an interval \eqn{\mathcal{Q} \in (0,1)}.
  
  This function uses numerical methods to solve for the
  end points of a HDR for a Beta density with user-specified shape
  parameters, via repeated calls to the functions \code{\link{dbeta}},
  \code{\link{pbeta}} and \code{\link{qbeta}}. The function
  \code{\link{optimize}} is used to find points \eqn{v} and \eqn{w}
  such that \deqn{f(v) = f(w)} subject to the constraint
  \deqn{\int_v^w f(\theta; \alpha, \beta) d\theta = p,}
  where \eqn{f(\theta; \alpha, \beta)} is a Beta density with shape
  parameters \eqn{\alpha} and \eqn{\beta}. 

  In the special case of \eqn{\alpha = \beta > 1}, the end points
  of a HDR  with content \eqn{p} are given by the \eqn{(1 \pm p)/2}
  quantiles of the Beta density, and are computed with the
  \code{\link{qbeta}} function. 

  Again note that the function will only compute a HDR for a unimodal
  Beta density, and exit with an error if \code{alpha<=1 | beta <=1}.
  Note that the uniform density results with \eqn{\alpha = \beta = 1},
  which does not have a unique HDR with content \eqn{0 < p <
  1}.  With shape parameters \eqn{\alpha<1} and \eqn{\beta>1} (or
  vice-versa, respectively), the Beta density is infinite at 0 (or 1,
  respectively), but still integrates to one, and so a HDR is still
  well-defined (but not implemented here, at least not yet).
  Similarly, with \eqn{0 < \alpha, \beta < 1} the Beta density is
  infinite at both 0 and 1, but integrates to one, and again a HDR of
  content \eqn{p<1} is well-defined in this case, but will be a set of
  two disjoint intervals (again, at present, this function does not
  cover this case).
}

\value{
  If the numerical optimization is successful an vector of length 2,
  containing \eqn{v} and \eqn{w}, defined above.    If the optimization
  fails for whatever reason, a vector of \code{NAs} is returned.
  
  The function will also produce a plot of the density with area under
  the density supported by the HDR shaded, if the user calls the
  function with \code{plot=TRUE}; the plot will appear on the current
  graphics device.

  Debugging messages are printed to the console if the \code{debug}
  logical flag is set to \code{TRUE}.  
}
\author{Simon Jackman \email{jackman@stanford.edu}.  Thanks to John
  Bullock who discovered a bug in an earlier version.}
\seealso{\code{\link{pbeta}}, \code{\link{qbeta}},
  \code{\link{dbeta}}, \code{\link{uniroot}}}
\examples{
betaHPD(4,5)
betaHPD(2,120)
betaHPD(120,45,p=.75,xlim=c(0,1))
}
\keyword{distribution}% at least one, from doc/KEYWORDS
