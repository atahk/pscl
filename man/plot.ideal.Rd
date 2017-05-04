\name{plot.ideal}
\alias{plot.ideal}
\alias{plot1d}
\alias{plot2d}
\title{plots an ideal object}
\description{
  Plot of the results of an ideal point estimation contained
  in an object of class \code{ideal}.
}

\usage{
\method{plot}{ideal}(x, conf.int=0.95, burnin=NULL, ...)

plot1d(x, d=1, conf.int=0.95, burnin=NULL,
       showAllNames = FALSE, ...)

plot2d(x, d1=1, d2=2, burnin=NULL,
       overlayCuttingPlanes=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{\link{ideal}}}
  \item{conf.int}{for "ideal" objects with 1 dimension estimated, the
    level of the confidence interval to plot around the
    posterior mean for each legislator.  If 2 or more dimensions were
    estimated, \code{conf.int} is ignored.}
  
  \item{d}{integer, which dimension to display in a 1d plot, if the
    object is a multidimensional ideal object}.
  
    \item{burnin}{of the recorded MCMC samples, how many to discard as
    burnin?  Default is \code{NULL}, in which case the value of
    \code{burnin} in the \code{\link{ideal}} object is used.}

  \item{showAllNames}{\code{\link{logical}}, if \code{TRUE}, the
    vertical axis will the names of all legislators. Default is
    \code{FALSE} to reduce clutter on typical-sized graph.}
  
  \item{d1}{integer, the number of the first dimension to plot when
    plotting multi-dimensional \code{\link{ideal}} objects.  This
    dimension will appear on the horizontal (x) axis.}
  
  \item{d2}{integer, the number of the second dimension to plot when
    plotting multi-dimensional \code{\link{ideal}} objects.  This dimension
    will appear on the vertical (y) axis.}
  
  \item{overlayCuttingPlanes}{logical, if \code{TRUE}, overlay the
    estimated bill-specific cutting planes}
  
  \item{...}{other parameters to be passed through to plotting functions.}
}


\details{ If the \code{ideal} object comes from fitting a \code{d=1}
  dimensional model, then \code{plot.ideal} plots the mean of the
  posterior density over each legislator's ideal point, accompanied by a
  \code{conf.int} confidence interval.  In this case, \code{plot.ideal}
  is simply a wrapper function to \code{\link{plot1d}}.

  If the \code{ideal} object has \code{d=2} dimensions, then
  \code{\link{plot2d}} is called, which plots the (estimated) mean of
  the posterior density of each legislator's ideal point (i.e., the
  ideal point/latent trait is a point in 2-dimensional Euclidean space,
  and the posterior density for each ideal point is a bivariate
  density).  Single dimension summaries of the estimated ideal points
  (latent traits) can be obtained for multidimensional
  \code{\link{ideal}} objects by passing the \code{\link{ideal}} object
  directly to \code{plot1d} with \code{d} set appropriately.

  If the \code{\link{ideal}} object has \code{d>2} dimensions, a
  scatterplot matrix is produced via \code{\link{pairs}}, with the
  posterior means of the ideal points (latent traits) plotted against
  one another, dimension by dimension.
  
  For unidimensional and two-dimensional models, if party information is
  available in the \code{rollcall} object contained in the \code{ideal}
  object, legislators from different parties are plotted in different
  colors.  If the \code{ideal} object has more than 2 dimensions,
  \code{plot.ideal()} produces a matrix of plots of the mean ideal
  points of each dimension against the posterior mean ideal points of
  the other dimensions.

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

\seealso{\code{\link{ideal}};  
  \code{\link{tracex}} for trace plots, a graphical aid useful in
diagnosing convergence of the MCMC algorithms.} 

\examples{
\dontrun{
data(s109)
id1 <- ideal(s109,
             d=1,
             normalize=TRUE,
             store.item=TRUE,
             maxiter=500,   ## short run for examples
             burnin=100,
             thin=10)  

plot(id1)


id2 <- ideal(s109,
             d=2,
             store.item=TRUE,
             maxiter=11e2,
             burnin=1e2,
             verbose=TRUE,
             thin=25)

plot(id2,overlayCuttingPlanes=TRUE)

id2pp <- postProcess(id2,
                     constraints=list(BOXER=c(-1,0),
                       INHOFE=c(1,0),
                       CHAFEE=c(0,.25)))

plot(id2pp,overlayCuttingPlanes=TRUE)
}
} 
\keyword{hplot}

