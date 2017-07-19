\name{summary.ideal}
\alias{summary.ideal}
\title{summary of an ideal object}

\description{
  Provides a summary of the output from ideal point estimation contained
  in an object of class \code{ideal}.
}

\usage{
\method{summary}{ideal}(object, prob=.95,
              burnin=NULL,
              sort=TRUE,
              include.beta=FALSE,...)
}

\arguments{

  \item{object}{an object of class \code{\link{ideal}}.}

  \item{prob}{scalar, a proportion between 0 and 1, the content of the highest posterior density (HPD) interval to compute for the parameters}
  
  \item{burnin}{of the recorded MCMC samples, how many to discard as
    burnin?  Default is \code{NULL}, in which case the value of
    \code{burnin} in the \code{\link{ideal}} object is used.}
  
  \item{sort}{logical, default is \code{TRUE}, indicating that the
    summary of the ideal points be sorted by the estimated posterior means
    (lowest to highest)}
  
  \item{include.beta}{whether or not to calculate summary statistics of
    beta, if beta is available.  If the item parameters were not stored
    in the \code{ideal} object, then \code{include.beta} is ignored.}
  
  \item{...}{further arguments passed to or from other functions}
  
}

\value{
  An item of class \code{summary.ideal} with elements:
  \item{object}{the name of the ideal object as an
    \code{\link[=eval]{unevaluated}} \code{\link{expression}}, produced
    by \code{match.call()$object}}
  \item{xm}{\code{n} by \code{d} matrix of posterior means for the ideal points}
  \item{xsd}{\code{n} by \code{d} matrix of posterior means for the ideal points}
  \item{xHDR}{\code{n} by 2 by \code{d} array of HDRs for the ideal points}  
  \item{bm}{\code{m} by \code{d+1} matrix of posterior means for the
      item parameters}
  \item{bsd}{\code{m} by \code{d+1} matrix of posterior standard deviation for the
      item parameters}
  \item{bHDR}{\code{m} by 2 by \code{d+1} array of HDRs for the item parameters}  
  \item{bSig}{a \code{\link{list}} of length \code{d}, each component a
    vector of length \code{m}, of mode \code{logical}, equal to
    \code{TRUE} if the corresponding discrimination parameter is
    distinguishible from zero; see Details.  If \code{store.item} was
    set to \code{FALSE} when \code{\link{ideal}} was invoked, then
    \code{bSig} is a list of length zero.}   
  \item{party.quant}{if party information is available through the
    \code{rollcall} object that was used to run \code{ideal}, then
    \code{party.quant} gives the posterior mean of the legislators'
    ideal points by party, by dimension.  If no party information is
    available, then \code{party.quant=NULL}.}
}

\details{The test of whether a given discrimination parameter is
  distinguishible from zero first checks to see if the two most extreme
  \code{quantiles} are symmetric around .5 (e.g., as are the default
  value of .025 and .975).  If so, the corresponding quantiles of the
  MCMC samples for each discrimination parameter are inspected to see if
  they have the same sign. If they do, then the corresponding
  discrimination parameter is flagged as distinguishible from zero;
  otherwise not.
}

\note{When specifying a value of \code{burnin} different from that used
  in fitting the \code{\link{ideal}} object, note a distinction
  between the iteration numbers of the stored iterations, and the
  number of stored iterations.  That is, the \code{n}-th iteration
  stored in an \code{\link{ideal}} object will not be iteration
  \code{n} if the user specified \code{thin>1} in the call to
  \code{\link{ideal}}.  Here, iterations are tagged with their
  iteration number.  Thus, if the user called \code{\link{ideal}} with
  \code{thin=10} and \code{burnin=100} then the stored iterations are
  numbered \code{100, 110, 120, ...}.  Any future subsetting via a
  \code{burnin} refers to this iteration number.}

\seealso{\code{\link{ideal}}}

\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}

\examples{
f <- system.file("extdata","id1.rda",package="pscl")
load(f)
summary(id1)

\dontrun{
data(s109)
cl2 <- constrain.legis(s109,
                       x=list("KENNEDY (D MA)"=c(-1,0),
                         "ENZI (R WY)"=c(1,0),
                         "CHAFEE (R RI)"=c(0,-.5)),
                       d=2)
id2Constrained <- ideal(s109,
                        d=2,
                        priors=cl2,      ## priors (w constraints)
                        startvals=cl2,   ## start value (w constraints)
                        store.item=TRUE,
                        maxiter=5000,
                        burnin=500,
                        thin=25)

summary(id2Constrained,
        include.items=TRUE)
}
}

\keyword{classes}

