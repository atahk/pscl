\name{constrain.items}
\alias{constrain.items}
\title{constrain item parameters in analysis of roll call data}

\description{
  Sets constraints on specified item parameters in Bayesian analysis of
  roll call data by generating appropriate priors and start values for
  Markov chain Monte Carlo iterations.
}

\usage{
constrain.items(obj, dropList = list(codes = "notInLegis", lop = 0),
                x, d = 1)
}

\arguments{
  \item{obj}{an object of class \code{\link{rollcall}}.}
  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}})
    indicating which voting decisions, legislators and/or roll calls are
    to be excluded from the subsequent analysis; see
    \code{\link{dropRollCall}} for details.}
  \item{x}{a \code{\link{list}} containing elements with names matching
    votes found in \code{dimnames(object$votes)[[2]]} (but after any
    subsetting specified by \code{dropList}).  Each component of the
    list must be a vector containing \code{d} elements, specifying the
    value to which the item discrimination parameters should be
    constrained, in each of the \code{d} dimensions.  The intercept or
    item difficultly parameter will not be constrained.}
  \item{d}{numeric, positive integer, the number of dimensions for which
    to set up the priors and start values.}
}

\details{\code{constrain.items} and its cousin,
  \code{\link{constrain.legis}} are usefully thought of as
  \dQuote{pre-processor} functions, generating priors \emph{and} start
  values for both the item parameters and the ideal points.  For the
  items specified in \code{x}, the prior mean for each dimension is set
  to the value given in \code{x}, and the prior precision for each
  dimension is set to \code{1e12} (i.e., a near-degenerate
  \dQuote{spike} prior).  For the other items, the priors are set to a
  mean of 0 and precision 0.01.  All of the ideal points are given
  normal priors with mean 0, precision 1.

  Start values are also generated for both ideal points and item
  parameters.  The start values for the items specified in \code{x} are
  set to the values specified in \code{x}.  The list resulting from
  \code{constrain.items} can then be given as the value for the
  parameters \code{priors} and \code{startvals} when \code{\link{ideal}}
  is run. The user is responsible for ensuring that a sufficient number
  of items are constrained such that when \code{\link{ideal}} is run,
  the model parameters are identified.

  \code{\link{dropRollCall}} is first called to generate the desired
  roll call matrix.  The entries of the roll call matrix are mapped to
  \code{c(0,1,NA)} using the \code{codes} component of the
  \code{\link{rollcall}} \code{object}.  See the discussion in the
  documentation of \code{\link{ideal}} for details on the generation of
  start values.
}

\value{
  a list with elements:
  \item{xp}{prior means for ideal points.  A matrix of dimensions number
    of legislators in \code{obj} by \code{d}.}
  \item{xpv}{prior meansprecisions for ideal points.  A matrix of dimensions number
    of legislators in \code{obj} by \code{d}.}
  \item{bp}{prior means for item parameters.  A matrix of dimensions number
    of items or votes in \code{obj} by \code{d+1}.}
  \item{bpv}{prior meansprecisions for item parameters.  A matrix of dimensions number
    of items or votes in \code{obj} by \code{d+1}.}
  \item{xstart}{start values for ideal points.  A matrix of dimensions number
    of legislators in \code{obj} by \code{d}.}
  \item{bstart}{start values for ideal points.  A matrix of dimensions number
    of items or votes in \code{obj} by \code{d+1}.}
}

\seealso{
  \code{\link{rollcall}}, \code{\link{ideal}},
  \code{\link{constrain.legis}}
}
\examples{
\dontrun{
data(s109)
f <- system.file("extdata","id1.rda",package="pscl")
load(f)
id1sum <- summary(id1,include.beta=TRUE)
suspect1 <- id1sum$bSig[[1]]=="95% CI overlaps 0"
close60 <- id1sum$bResults[[1]][,"Yea"] < 60
close40 <- id1sum$bResults[[1]][,"Yea"] > 40
suspect <- suspect1 & close60 & close40
id1sum$bResults[[1]][suspect,]
suspectVotes <- dimnames(id1sum$bResults[[1]][suspect,])[[1]]


## constraints on 2d model,
## close rollcall poorly fit by 1d model
## serves as reference item for 2nd dimension

cl <- constrain.items(s109,
                      x=list("2-150"=c(0,7),
                        "2-169"=c(7,0)),
                      d=2)

id1Constrained <- ideal(s109,
                        d=2,
                        meanzero=TRUE,
                        priors=cl,
                        startvals=cl,
                        maxiter=1e5,
                        burnin=1e3,
                        thin=1e2)
summary(id1Constrained,include.beta=TRUE)
}
}
\keyword{datagen}

