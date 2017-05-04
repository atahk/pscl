\name{constrain.legis}
\alias{constrain.legis}
\title{constrain legislators' ideal points in analysis of roll call data}

\description{
  Sets constraints on specified legislators for ideal
  point estimation by generating appropriate priors and start values.
}

\usage{
constrain.legis(obj, dropList = list(codes = "notInLegis", lop = 0),
                 x, d = 1)

}

\arguments{
  \item{obj}{an object of class \code{\link{rollcall}}.}

  \item{dropList}{a \code{\link{list}} (or \code{\link{alist}})
    indicating which voting decisions, legislators and/or roll calls are
    to be excluded from the subsequent analysis; see
    \code{\link{dropRollCall}} for details.}

  \item{x}{a \code{\link{list}} containing elements with names partially
    matching legislators found in \code{dimnames(object$votes)[[1]]} (but after
    any sub-setting specified by \code{dropList}).  Each element must be a
    vector containing \code{d} elements, specifying the value to which
    the ideal point should be constrained in each of \code{d}
    dimensions. \code{x} must have at least \code{d+1} components; i.e.,
    supplying a necessary (but not sufficient) set of constraints for
    global identification of the parameters of a \code{d}-dimensional
    item-response model, see Details.}

  \item{d}{the number of dimensions for which to set up the priors and
    start values.}
}

\details{
  \code{constrain.items} and its cousin,
  \code{\link{constrain.legis}} are usefully thought of as
  \dQuote{pre-processor} functions, implementing identification
  constraints for the ideal point model by generating priors \emph{and} start
  values for both the item parameters and the ideal points.

  For the legislators specified in \code{x}, the prior mean for each
  dimension is set to the specified value and the prior precision for
  each dimension is set to \code{1e12} (i.e., a near-degenerate
  \dQuote{spike} prior, and, for all practical purposes, constraining
  that parameter to a fixed value).  For the other legislators, the
  priors on their ideal points are set to a mean of 0 and a small
  precision of .01, corresponding to a prior variance of 100, or a prior
  95 percent confidence interval of -20 to 20.  All of the item
  parameter priors are set to mean 0, precision 0.01.

  Start values are also generated for both ideal points and item
  parameters.  The start values for the legislators named in
  \code{x} are set to the values specified in \code{x}.  The list
  resulting from \code{constrain.legis} can then be given as the value
  for the parameters \code{priors} and \code{startvals} when
  \code{\link{ideal}} is run. \code{constrain.legis} requires that
  \code{d+1} constraints be specified; if the constrained ideal points
  points are linearly independent, then the parameters of the
  item-response model are (at least locally) identified.  For instance,
  when fitting a 1 dimensional model, constraining the ideal points of
  two legislators is sufficient to globally identify the model
  parameters.

 \code{\link{dropRollCall}} is first called to generate the desired
  roll call matrix.  The entries of the roll call matrix are mapped to
  \code{c(0,1,NA)} using the \code{codes} component of the
  \code{\link{rollcall}} \code{object}.  See the discussion in the
  documentation of \code{\link{ideal}} for details on the generation of
  start values.

}


\value{ a list with elements:

  \item{xp}{prior means for ideal points.  A matrix of dimensions number
    of legislators in \code{rc} by \code{d}.}
  
  \item{xpv}{prior meansprecisions for ideal points.  A matrix of
    dimensions number of legislators in \code{rc} by \code{d}.}
  
  \item{bp}{prior means for item parameters.  A matrix of
    dimensions number of items or votes in \code{rc} by \code{d+1}.}
  
  \item{bpv}{prior meansprecisions for item parameters.  A matrix of
    dimensions number of items or votes in \code{rc} by \code{d+1}.}
  
  \item{x}{start values for ideal points.  A matrix of dimensions
    number of legislators in \code{rc} by \code{d}.}
  
  \item{b}{start values for ideal points.  A matrix of dimensions
    number of items or votes in \code{rc} by \code{d+1}.}

}

\seealso{
  \code{\link{rollcall}}, \code{\link{ideal}},
  \code{\link{constrain.items}}.
  See \code{\link{pmatch}} on how supplied names are matched against the
  names in the \code{\link{rollcall}} object.
} 


\examples{
data(s109)
cl <- constrain.legis(s109,
                      x=list("KENNEDY"=-1,
                        "ENZI"=1),
                      d=1)

\dontrun{
## too long for examples
id1Constrained <- ideal(s109,
                       d=1,
                       priors=cl,      ## use cl
                       startvals=cl,   ## use cl
                       maxiter=5000,
                       burnin=500,
                       thin=25)
summary(id1Constrained)

cl2 <- constrain.legis(s109,
                       x=list("KENNEDY"=c(-1,0),
                         "ENZI"=c(1,0),
                         "CHAFEE"=c(0,-.5)),
                       d=2)


id2Constrained <- ideal(s109,
                        d=2,
                        priors=cl2,      ## priors (w constraints)
                        startvals=cl2,   ## start value (w constraints)
                        store.item=TRUE,
                        maxiter=5000,
                        burnin=500,
                        thin=25)
summary(id2Constrained,include.items=TRUE)
}
}
\keyword{datagen}

