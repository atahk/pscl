\name{postProcess}
\alias{postProcess}
\title{remap MCMC output via affine transformations}

\description{
  Remap the MCMC iterates in an \code{\link{ideal}} object via an affine
  transformation, imposing identifying restrictions ex post (aka
  post-processing).} 

\usage{
postProcess(object, constraints="normalize", debug = FALSE)
}

\arguments{ \item{object}{an object of class \code{\link{ideal}}}
  \item{constraints}{list of length \code{d+1}, each component providing
  a set of \code{d} restrictions, where \code{d} is the dimension of the
  fitted \code{\link{ideal}} model; or the character string
  \code{normalize} (default).  If a list, the name of each component should
  uniquely match a legislator/subject's name.  See Details.}
  \item{debug}{logical flag for verbose output, used for debugging} }

\details{
  Item-response models are unidentified without restrictions on the
  underlying parameters.  Consider the \code{d=1} dimensional case. The
  model is
  \deqn{P(y_{ij} = 1) = F(x_i \beta_j - \alpha_j)}{%
    Pr(y_[ij] = 1) = F(x_i b_j - a_j).
  }
  Any
  linear transformation of the latent traits, say,
  \deqn{x^* = mx + c}{%
    x* = mx + c
  }
  can be exactly offset by applying the appropriate linear
  transformations to the item/bill parameters, meaning that there is no
  unique set of values for the model parameters that will maximize the
  likelihood function.  In higher dimensions,
  the latent traits can also be transformed via any arbitrary rotation,
  dilation and translation, with offsetting transformations applied to
  the item/bill parameters.

  One strategy in MCMC is to ignore the lack of identification at run
  time, but apply identifying restrictions ex post,
  \dQuote{post-processing} the MCMC output, iteration-by-iteration.  In
  a \code{d}-dimensional IRT model, a sufficient condition for global
  identification is to fix \code{d+1} latent traits, provided the
  constrained latent traits span the \code{d} dimensional latent space.
  This function implements this strategy.  The user supplies a set of
  constrained ideal points in the \code{constraints} list.  The function
  then processes the MCMC output in the \code{\link{ideal}}
  \code{object}, finding the transformation that maps the current
  iteration's sampled values for \code{x} (latent traits/ideal points)
  into the sub-space of identified parameters defined by the fixed
  points in \code{constraints}; i.e., what is the affine transformation
  that maps the unconstrained ideal points into the constraints?  Aside
  from minuscule numerical inaccuracies resulting from matrix inversion
  etc, this transformation is exact: after post-processing, the
  \code{d+1} constrained points do not vary over the MCMC iterations.
  The remaining \code{n-d-1} ideal points are subject to (posterior)
  uncertainty; the \dQuote{random tour} of the joint parameter space of
  these parameters produced by the MCMC algorithm has been mapped into a
  subspace in which the parameters are globally identified.

  If the \code{\link{ideal}} object was produced with \code{store.item}
  set to \code{TRUE}, then the item parameters are also post-processed,
  applying the inverse transformation.  Specifically, recall that the
  IRT model is \deqn{P(y_{ij} = 1) = F(x_i'\beta_j)}{%
    Pr(y_[ij] = 1) = F(x_i' b_j)}
  where in this formulation \eqn{x_i} is a vector of
  length \code{d+1}, including a \code{-1} to put a constant term into
  the model (i.e., the intercept or difficulty parameter
  is part of \eqn{\beta_j}{beta_j}).  Let \eqn{A} denote the
  non-singular, \code{d+1}-by-\code{d+1} matrix that maps the \eqn{x}
  into the space of identified parameters.  Recall that this
  transformation is computed iteration by iteration.  Then each
  \eqn{x_i} is transformed to \eqn{x^*_i = Ax_i}{x*_i = Ax_i} and
  \eqn{\beta_j}{b_j} is transformed to \eqn{\beta_j^* = A^{-1}
    \beta_j}{b_j^* = A^(-1) b_j}, \eqn{i = 1, \ldots, n; j = 1,
    \ldots, m}.
  
  Local identification can be obtained for a one-dimensional model by
  simply imposing a normalizing restriction on the ideal points: this
  normalization (mean zero, standard deviation one) is the default
  behavior, but (a) is only sufficient for local identification when
  the \code{rollcall} object was fit with \code{d=1}; (b) is not
  sufficient for even local identification when \code{d>1}, with
  further restrictions required so as to rule out other forms of
  invariance (e.g., translation, or "dimension-switching", a
  phenomenon akin to label-switching in mixture modeling).

  The default is to impose dimension-by-dimension normalization with
  respect to the means of the marginal posterior densities of the
  ideal points, such that the these means (the usual Bayes estimates
  of the ideal points) have mean zero and standard deviation one
  across legislators.  An offsetting transformation is applied to the
  items parameters as well, if they are saved in the \code{ideal} object.  

  Specifically, in one-dimension, the
  two-parameter IRT model is \deqn{P(y_{ij} = 1) = F(x_i \beta_j -
  \alpha_j).}{Pr(y_[ij] = 1) = F(x_i b_j - a_j).}  If we normalize the
  \eqn{x_i} to \eqn{x*_i = (x_i - c)/m} then the offsetting
  transformations for the item/bill parameters are \eqn{\beta_j^* =
  \beta_j m}{b*_j = b_j m} and \eqn{\alpha_j^* = \alpha_j -
  c\beta_j}{a*_j = a_j - cb_j}.
}

\value{
  An object of class \code{\link{ideal}}, with components suitably
  transformed and recomputed (i.e., \code{x} is transformed and
  \code{xbar} recomputed, and if the
  \code{\link{ideal}} object was fit with \code{store.item=TRUE},
  \code{beta} is transformed and \code{betabar} is recomputed).
}

\references{
Hoff, Peter, Adrian E. Raftery and Mark S. Handcock. 2002.
Latent Space Approaches to Social Network Analysis. \emph{Journal
  of the American Statistical Association} 97:1090--1098.

Edwards, Yancy D. and Greg M. Allenby. 2003.
Multivariate Analysis of Multiple Response Data. \emph{Journal of
  Marketing Research} 40:321--334.

Rivers, Douglas. 2003.  \dQuote{Identification of Multidimensional
  Item-Response Models.} Typescript.  Department of Political Science,
Stanford University.
}

\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}

\note{Applying transformations to obtain identification can
  sometimes lead to surprising results.  Each data point makes the same
  likelihood contributions with either the identified or unidentified
  parameters.  But, in general, predictions generated with the
  parameters set to their posterior means will differ depending on
  whether one uses the identified subset of parameters or the
  unidentified parameters.  For this reason, caution should be used when
  using a function such as \code{\link{predict}} after post-processing
  output from \code{\link{ideal}}.  A better strategy is to compute the
  estimand of interest at each iteration and then take averages over
  iterations.

  When specifying a value of \code{burnin} different from that used in
  fitting the \code{\link{ideal}} object, note a distinction between
  the iteration numbers of the stored iterations, and the number of
  stored iterations.  That is, the \code{n}-th iteration stored in an
  \code{\link{ideal}} object will not be iteration \code{n} if the
  user specified \code{thin>1} in the call to \code{\link{ideal}}.
  Here, iterations are tagged with their iteration number.  Thus, if
  the user called \code{\link{ideal}} with \code{thin=10} and
  \code{burnin=100} then the stored iterations are numbered \code{100,
    110, 120, ...}.  Any future subsetting via a \code{burnin} refers to
  this iteration number.
}
  
\examples{
data(s109)
f = system.file("extdata",package="pscl","id1.rda")
load(f)

id1Local <- postProcess(id1)    ## default is to normalize
summary(id1Local)

id1pp <- postProcess(id1,
                     constraints=list(BOXER=-1,INHOFE=1))
summary(id1pp)

## two-dimensional fit
f = system.file("extdata",package="pscl","id2.rda")
load(f)

id2pp <- postProcess(id2,
                     constraints=list(BOXER=c(-1,0),
                       INHOFE=c(1,0),
                       CHAFEE=c(0,.25)))

tracex(id2pp,d=1:2,
       legis=c("BOXER","INHOFE","COLLINS","FEINGOLD","COLEMAN",
         "CHAFEE","MCCAIN","KYL"))
}
\keyword{models}
