\name{tracex}
\alias{tracex}
\title{trace plot of MCMC iterates, posterior density of legislators'
  ideal points} 

\description{
  Produces a trace plot of the MCMC samples from the posterior density of
  legislators' \code{\link{ideal}} points.
}

\usage{
tracex(object, legis=NULL, d=1, conf.int=0.95,
       multi = FALSE, burnin=NULL,span=.25,
       legendLoc="topright")
}

\arguments{
  \item{object}{an object of class \code{ideal}.}
  \item{legis}{a vector of either the names of legislators (or
  \code{\link[=pmatch]{partial matches}} of the names as given in the
  \code{\link{dimnames}} of \code{object$x}.}
  \item{d}{numeric, either a scalar or a vector of length two, the
  dimension(s) to be traced.} 
  \item{conf.int}{numeric, the level of the confidence interval on the
  posterior mean to be plotted.}
  \item{multi}{logical, multiple plotting panels, one per legislators? If \code{FALSE} (default) and \code{length(d)==2},
  display traces for all selected legislators' ideal points on
  the one plot.}
 \item{burnin}{of the recorded MCMC samples, how many to discard as
    burnin?  Default is \code{NULL}, in which case the value of
    \code{burnin} in the \code{\link{ideal}} object is used.}
  \item{span}{numeric, a proportion, the \code{span} to be used when
    calling \code{loess} to generate a moving average for trace plots when
    \code{d=1}}
  \item{legendLoc}{numeric or character, and possibly a vector,
  specifying where to place the legend when \code{d=1}; setting
  \code{legendLoc=NULL} will suppress the legend for all requested trace
  plots}
}

\details{Produces a trace plot showing the history of the MCMC
  iterations for the ideal point of each of the legislators (partially)
  named in \code{legis}.  For \code{d=1}, each trace plot includes a
  trace over iterations, the cumulative mean, a moving average, the
  MCMC-based estimate of the mean of the posterior, and a confidence
  interval (specified by \code{conf.int}) around the mean of the
  posterior (using the estimated \code{\link[=quantile]{quantiles}}) of
  the respective MCMC iterates).  All of these values are calculated
  discarding the initial \code{burnin} iterations.

  When \code{d} is a vector of length two, a 2-dimensional trace plot is
  displayed, with the \code{d[1]} dimension on the horizontal axis, and
  the \code{d[2]} dimension on the vertical axis.

  When \code{d=1}, a legend will be placed on the plot; the option
  \code{legendLoc} controls the placing of the legend.  \code{legendLoc}
  may be a vector, specifying a unique legend location for each
  requested trace plots.  If
  \code{legendLoc} is of length 1, it will be \code{\link{rep}}licated
  to have length equal to the number of requested trace plots.
}



\seealso{\code{\link{ideal}}; \code{\link{pmatch}} for matching
  legislators' names.  See \code{\link{legend}} for valid options to
  \code{legendLoc}.} 

\examples{
data(s109)
f <- system.file("extdata","id1.rda",package="pscl")
load(f)
tracex(id1,legis="KENN")

## n.b., no such legislator named Thomas Bayes
tracex(id1,legis=c("KENN","BOX","KYL","Thomas Bayes"))


f <- system.file("extdata","id2.rda",package="pscl")
load(f)

tracex(id2,d=1,legis=c("KENNEDY","BOXER","KYL","Thomas Bayes"))
tracex(id2,d=2,legis=c("KENNEDY","BOXER","KYL","Thomas Bayes"))
tracex(id2,d=1:2,
       legis=c("KENNEDY","BOXER","KYL","Thomas Bayes"))

## partial matching
tracex(id2,d=1:2,
       legis=c("KENN","BOX","BID","SNO","SPEC","MCCA","KYL",
               "Thomas Bayes"),
       multi=TRUE)
}

\keyword{hplot}

