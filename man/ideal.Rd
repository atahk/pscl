\name{ideal}
\alias{ideal}
\title{analysis of educational testing data and roll call data with IRT models, via Markov
  chain Monte Carlo methods}

\description{
  Analysis of \code{rollcall} data via the spatial voting model;
  equivalent to a 2 parameter item-response model to educational testing data.  Model fitting via Markov chain Monte Carlo (MCMC).
}

\usage{
ideal(object, codes = object$codes,
      dropList = list(codes = "notInLegis", lop = 0),
      d = 1, maxiter = 10000, thin = 100, burnin = 5000,
      impute = FALSE,
      normalize = FALSE,
      meanzero = normalize,
      priors = NULL, startvals = "eigen",
      store.item = FALSE, file = NULL,
      verbose=FALSE, use.voter=NULL)
}

\arguments{
  \item{object}{an object of class \code{\link{rollcall}}}
  
  \item{codes}{a \code{\link{list}} describing the types of voting
    decisions in the roll call matrix (the \code{votes} component of the
    \code{\link{rollcall}} \code{object}); defaults to 
    \code{object$codes}, the codes in the rollcall object.} 
  
  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}})
    listing voting decisions, legislators and/or votes to be dropped
    from the analysis; see \code{\link{dropRollCall}} for details.}
  
  \item{d}{numeric, (small) positive integer (default = 1), dimensionality of the ability space (or "policy space" in the rollcall context).}

  \item{maxiter}{numeric, positive integer, multiple of \code{thin}, number of MCMC iterations}

  \item{thin}{numeric, positive integer, thinning interval used for
    recording MCMC iterations.}

  \item{burnin}{number of MCMC iterations to run before recording.  The
    iteration numbered \code{burnin} will be recorded.  Must be a
    multiple of \code{thin}.}

  \item{impute}{\code{\link{logical}}, whether to treat missing entries
    of the rollcall matrix as missing at random, sampling from the
    predictive density of the missing entries at each MCMC iteration.}
  	    
  \item{normalize}{\code{\link{logical}}, impose identification with
    the constraint that the ideal points have mean zero and
    standard deviation one, in each dimension.  For one dimensional models this option is sufficient to
    locally identify the model parameters.
    See Details.}

  \item{meanzero}{to be deprecated/ignored; use \code{normalize} instead.}
  
  \item{priors}{a \code{list} of parameters (means and variances)
    specifying normal priors for the legislators' ideal points. The
    default is \code{NULL}, in which case the normal priors used have mean zero and
    precision 1 for the ideal points (ability parameters) and mean zero and
    precision .04 (variance 25) for the bill parameters (item discrimination and difficulty parameters).  If not \code{NULL}, \code{priors} must be a
    \code{list} with as many as four named components \code{xp, xpv, bp,
      bpv}:
    \describe{
    \item{\code{xp}}{a \code{n} by \code{d} matrix
      of prior \emph{means} for the legislators' ideal points;
      or alternatively, a scalar, which will be replicated to fill a \code{n} by \code{d} matrix.}
    \item{\code{xpv}}{a \code{n} by \code{d} matrix of prior
      \emph{precisions} (inverse variances);
      or alternatively, a scalar, which will be replicated to fill a \code{n}
      by \code{d} matrix.} 
    \item{\code{bp}}{a \code{m} by \code{d+1} matrix of prior means for the
      item parameters (with the item difficulty parameter coming last);
      or alternatively, a scalar, which will be replicated to fill a \code{m}
      by \code{d+1} matrix.}
    \item{\code{bpv}}{a \code{m} by \code{d+1} matrix of prior
      precisions for the item parameters;
      or alternatively, a scalar, which will be replicated to fill a \code{m}
      by \code{d+1} matrix.}
    }
    None of the components should contain \code{NA}.  If any
    of the four possible components are not provided, then the
    corresponding component of \code{priors} is assigned using the default
    values described above.}

  \item{startvals}{either a string naming a method for generating start
    values, valid options are \code{"eigen"} (the default),
    \code{"random"} or a \code{list} containing start values for
    legislators' ideal points and item parameters.  See Details.}

  \item{store.item}{\code{\link{logical}}, whether item discrimination
    parameters should be stored.  Storing item discrimination parameters
    can consume a large amount of memory.  These need to be stored for
    prediction; see \code{\link{predict.ideal}}.}

  \item{file}{string, file to write MCMC output.  Default is
    \code{NULL}, in which case MCMC output is stored in memory.  Note
    that post-estimation commands like \code{plot} will not work unless
    MCMC output is stored in memory.}

  \item{verbose}{logical, default
    is \code{FALSE}, which generates relatively little output to the R
    console during execution.}

  \item{use.voter}{A vector of logicals of length \code{n} controlling
    which legislators' vote data informs item parameter
    estimates. Legislators corresponding to \code{FALSE} entries will
    not have their voting data included in updates of the item
    parameters. The default value of \code{NULL} will run the standard
    ideal-point model, which uses all legislators in updating item
    parameters. See Jessee (2016).}
}

\details{The function fits a \code{d}+1 parameter item-response model to
  the roll call data object, so in one dimension the model reduces
  to the two-parameter item-response model popular in educational testing.
  See References.
  
  \strong{Identification}: The model parameters are \strong{not
    identified} without the user supplying some restrictions on the
    model parameters; i.e., translations, rotations and re-scalings of
    the ideal points are observationally equivalent, via offsetting
    transformations of the item parameters.  It is the user's
    responsibility to impose these identifying restrictions if
    desired. The following brief discussion provides some guidance.

  For one-dimensional models (i.e., \code{d=1}), a simple route to
  identification is the \code{normalize} option, by imposing the restriction that the means of the posterior densities of the ideal points (ability parameters) have mean zero and standard deviation one, across legislators (test-takers).  This normalization supplies
  \emph{local} identification (that is, identification up to a 180 degree rotation of
  the recovered dimension). 
  
  Near-degenerate \dQuote{spike} priors
  (priors with arbitrarily large precisions) or the
  \code{constrain.legis} option on any two legislators' ideal points
  ensures \emph{global} identification in one dimension.

  Identification in higher dimensions can be obtained by supplying
  fixed values for \code{d+1} legislators' ideal points, provided the
  supplied fixed points span a \code{d}-dimensional space (e.g., three
  supplied ideal points form a triangle in \code{d=2} dimensions), via
  the \code{\link{constrain.legis}} option. In this case the function
  defaults to vague normal priors on the unconstrained ideal points, but at each iteration the sampled
  ideal points are transformed back into the space of identified
  parameters, applying the linear transformation that maps the
  \code{d+1} fixed ideal points from their sampled values to their
  fixed values.  Alternatively, one can impose
  restrictions on the item parameters via
  \code{\link{constrain.items}}. See the examples in the documentation
  for the \code{\link{constrain.legis}} and
  \code{\link{constrain.items}}.

  Another route to identification is via \emph{post-processing}.  That
  is, the user can run \code{ideal} without any identification
  constraints. This does not pose any formal/technical problem in a
  Bayesian analysis.  The fact that the posterior density may have
  mulitple modes doesn't imply that the posterior is improper or that
  it can't be explored via MCMC methods. -- but then use the function
  \code{\link{postProcess}} to map the MCMC output from the space of
  unidentified parameters into the subspace of identified parameters.
  See the example in the documentation for the
  \code{\link{postProcess}} function.

  When the \code{normalize} option is set to \code{TRUE}, an
  unidentified model is run, and the \code{ideal} object is
  post-processed with the \code{normalize} option, and then returned
  to the user (but again, note that the \code{normalize} option is
  only implemented for unidimensional models).
  
  \strong{Start values}.  Start values can be supplied by the user, or
  generated by the function itself. 

  The default method, corresponding to \code{startvals="eigen"}, first
  forms a \code{n}-by-\code{n} correlation matrix from the
  double-centered roll call matrix (subtracting row means, and column
  means, adding in the grand mean), and then extracts the first
  \code{d} principal components (eigenvectors), scaling the
  eigenvectors by the square root of their corresponding eigenvector.
  If the user is imposing constraints on ideal points (via
  \code{\link{constrain.legis}}), these constraints are applied to the
  corresponding elements of the start values generated from the
  eigen-decomposition.  Then, to generate start 
  values for the rollcall/item parameters, a series of
  \code{\link[=family]{binomial}} \code{\link[=glm]{glms}} are
  estimated (with a probit \code{\link[=make.link]{link}}), one for
  each rollcall/item, \eqn{j = 1, \ldots, m}.  The votes on the
  \eqn{j}-th rollcall/item are binary responses (presumed to be
  conditionally independent given each legislator's latent
  preference), and the (constrained or unconstrained) start values for
  legislators are used as predictors. The estimated coefficients from
  these probit models are used as start values for the item
  discrimination and difficulty parameters (with the intercepts from
  the probit GLMs multiplied by -1 so as to make those coefficients
  difficulty parameters).

  The default \code{eigen} method generates extremely good start
  values for low-dimensional models fit to recent U.S. congresses,
  where high rates of party line voting result in excellent fits from
  low dimensional models.  The \code{eigen} method may be
  computationally expensive or lead to memory errors for
  \code{rollcall} objects with large numbers of legislators.

  The \code{random} method generates start values via iid sampling
  from a N(0,1) density, via \code{\link{rnorm}}, imposing any
  constraints that may have been supplied via
  \code{\link{constrain.legis}}, and then uses the probit method
  described above to get start values for the rollcall/item
  parameters.
  
  If \code{startvals} is a \code{list}, it must contain the named
  components \code{x} and/or \code{b}, or named components that
  (uniquely) begin with the letters \code{x} and/or \code{b}.  The
  component \code{x} must be a vector or a matrix of dimensions equal to
  the number of individuals (legislators) by \code{d}.  If supplied,
  \code{startvals$b} must be a matrix with dimension number of items
  (votes) by \code{d}+1.  The \code{x} and \code{b} components cannot
  contain \code{NA}.  If \code{x} is not supplied when \code{startvals}
  is a list, then start values are generated using the default
  \code{eiegn} method described above, and start values for the
  rollcall/item parameters are regenerated using the probit method,
  ignoring any user-supplied values in \code{startvals$b}.  That is,
  user-supplied values in \code{startvals$b} are only used when
  accompanied by a valid set of start values for the ideal points in
  \code{startvals$x}.

  \strong{Implementation via Data Augmentation}.  The MCMC algorithm
  for this problem consists of a Gibbs sampler for the ideal points
  (latent traits)
  and item parameters, conditional on latent data \eqn{y^*}, generated
  via a data augmentation (DA) step.  That is, following Albert (1992) and
  Albert and Chib (1993), if \eqn{y_{ij} = 1} we sample from the
  truncated normal density
  \deqn{y_{ij}^* \sim N(x_i' \beta_j - \alpha_j, 1)\mathcal{I}(y_{ij}^* \geq 0)}
  and for \eqn{y_{ij}=0} we sample 
  \deqn{y_{ij}^* \sim N(x_i' \beta_j - \alpha_j, 1)\mathcal{I}(y_{ij}^* < 0)}
  where \eqn{\mathcal{I}} is an indicator function evaluating to one
  if its argument is true and zero otherwise.  Given the
  latent \eqn{y^*}, the conditional distributions for \eqn{x} and
  \eqn{(\beta,\alpha)} are extremely simple to sample from; see the
  references for details.  

  This data-augmented Gibbs sampling strategy is easily implemented,
  but can sometimes require many thousands of samples in order to
  generate tolerable explorations of the posterior densities of the
  latent traits, particularly for legislators with short and/or
  extreme voting histories (the equivalent in the educational testing
  setting is a test-taker who gets almost every item right or wrong).

  % The MCMC algorithm can generate better performance
  % via a parameter expansion strategy usually referred to as \emph{marginal
  % data augmentation} (e.g., van Dyk and Meng 2001).  The idea is to
  % introduce a additional working parameter into the MCMC sampler that
  % has the effect of improving the performance of the sampler in the
  % sub-space of parameters of direct interest.  In this case we
  % introduce a variance parameter \eqn{\sigma^2} for the latent data;
  % in the DA algorithm of Albert and Chib (1993) --- and in any conventional
  % probit analysis --- this parameter is set
  % to 1.0 for identification.  In the MDA approach we carry this
  % (unidentified) parameter into the DA stage of the algorithm with an
  % improper prior, \eqn{p(\sigma^2) \propto \sigma^{-2}},
  % generating \eqn{y^*} that exhibit bigger moves from iteration to
  % iteration, such that in turn the MCMC algorithm displays better
  % mixing with respect to the identified parameters of direct interest,
  % \eqn{x} and \eqn{(\beta,\alpha)} than the mixing obtained from
  % the Gibbs-with-DA MCMC algorithm.   The MDA algorithm is the default
  % in \code{ideal}, but Gibbs-with-DA can be implemented by setting
  % \code{mda=FALSE} in the call to \code{ideal}.
}
  
\value{a \code{\link{list}} of class \code{ideal} with named components
  
  \item{n}{\code{\link{numeric}}, integer, number of legislators in the
    analysis, after any subseting via processing the \code{dropList}.}

  \item{m}{\code{\link{numeric}}, integer, number of rollcalls in roll
    call matrix, after any subseting via processing the \code{dropList}.}

  \item{d}{\code{\link{numeric}}, integer, number of dimensions
    fitted.}

  \item{x}{a three-dimensional \code{\link{array}} containing the MCMC
    output with respect to the 
    the ideal point of each legislator in each dimension.
    The three-dimensional array is in iteration-legislator-dimension
    order.  The iterations run from \code{burnin} to \code{maxiter}, at
    an interval of \code{thin}.}

   \item{beta}{a three-dimensional \code{\link{array}} containing the
    MCMC output for the item parameters.  The three-dimensional array
    is in iteration-rollcall-parameter order.  The iterations run from
    \code{burnin} to \code{maxiter}, at an interval of \code{thin}.
    Each rollcall has \code{d+1} parameters, with the
    item-discrimination parameters stored first, in the first \code{d}
    components of the 3rd dimension of the \code{beta} array; the
    item-difficulty parameter follows in the final \code{d+1}
    component of the 3rd dimension of the \code{beta} array.}

  \item{xbar}{a \code{n} by \code{d} \code{\link{matrix}} containing the means of the
    MCMC samples for the ideal point of each legislator in each dimension,
    using iterations \code{burnin} to \code{maxiter}, at an interval of
    \code{thin}.}
  
  \item{betabar}{a \code{m} by \code{d+1} \code{\link{matrix}} containing the means of
    the MCMC samples for the item-specific parameters, using iterations
    \code{burnin} to \code{maxiter}, at an interval of \code{thin}.}

	\item{args}{calling arguments, evaluated in the frame calling \code{ideal}.}

  \item{call}{an object of class \code{\link{call}}, containing
    the arguments passed to \code{ideal} as unevaluated expressions or values (for functions arguments that evaluate to scalar integer or logical such as \code{maxiter}, \code{burnin}, etc).}
}

\references{
  Albert, James. 1992. Bayesian Estimation of normal ogive item
  response curves using Gibbs sampling. \emph{Journal of Educational
  Statistics}. 17:251-269.

  Albert, James H. and Siddhartha Chib.  1993.  Bayesian Analysis of
  Binary and Polychotomous Response Data.  \emph{Journal of the
  American Statistical Association}. 88:669-679.

  Clinton, Joshua, Simon Jackman and Douglas Rivers. 2004. The
  Statistical Analysis of Roll Call Data.  \emph{American Political
  Science Review}.  98:335-370.

  Jackman, Simon. 2009.  \emph{Bayesian Analysis for the Social
  Sciences}.  Wiley: Hoboken, New Jersey.

  Jessee, Stephen. 2016. (How) Can We Estimate the Ideology
  of Citizens and Political Elites on the Same Scale?
  \emph{American Journal of Political Science}.

  Patz, Richard J. and Brian W. Junker. 1999.  A Straightforward
  Approach to Markov Chain Monte Carlo Methods for Item Response
  Models. \emph{Journal of Education and Behavioral
  Statistics}. 24:146-178.
  
  Rivers, Douglas. 2003.  \dQuote{Identification of Multidimensional
  Item-Response Models.} Typescript.  Department of Political Science,
  Stanford University.

  van Dyk, David A and Xiao-Li Meng.  2001.  The art of data
  augmentation (with discussion).  \emph{Journal of Computational and
  Graphical Statistics}. 10(1):1-111.
}

\author{Simon Jackman \email{jackman@stanford.edu}, with help from Christina
  Maimone and Alex Tahk.}

\seealso{
  \code{\link{rollcall}}, \code{\link{summary.ideal}},
  \code{\link{plot.ideal}}, \code{\link{predict.ideal}}.
  \code{\link{tracex}} for graphical display of MCMC iterative
  history.

  \code{\link{idealToMCMC}} converts the MCMC iterates in an
  \code{ideal} object to a form that can be used by the \code{coda} library.

  \code{\link{constrain.items}} and
  \code{\link{constrain.legis}} for implementing identifying
  restrictions.

  \code{\link{postProcess}} for imposing identifying restrictions
  \emph{ex post}.
  
  \code{\link[MCMCpack:MCMCirt1d]{MCMCirt1d}} and 
  \code{\link[MCMCpack:MCMCirtKd]{MCMCirtKd}} in the \pkg{MCMCpack}
  package provide similar functionality to \code{ideal}.
  
}

\examples{
\dontrun{
## long run, many iterations
data(s109)
n <- dim(s109$legis.data)[1]
x0 <- rep(0,n)
x0[s109$legis.data$party=="D"] <- -1
x0[s109$legis.data$party=="R"] <- 1

id1 <- ideal(s109,
             d=1,
             startvals=list(x=x0),
             normalize=TRUE,
             store.item=TRUE,
             maxiter=260E3,
             burnin=10E3,
             thin=100)  
}
}
\keyword{models}
