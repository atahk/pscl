\name{plot.seatsVotes}
\alias{plot.seatsVotes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot seats-votes curves}
\description{
  Plot seats-votes curves produced by \code{seatsVotes}
}
\usage{
\method{plot}{seatsVotes}(x, type = c("seatsVotes", "density"),
                legend = "bottomright", transform=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{\link{seatsVotes}}}
  \item{type}{character, partially matching the options above; see details}
  \item{legend}{where to put the legend when plotting with
    \code{type="seatsVotes"}}
  \item{transform}{logical, whether to transform the vote shares for
    \code{type="density"}; see Details} 
  \item{\dots}{arguments passed to or from other functions (e.g.,
    options for the \code{\link{density}} function, when
    \code{type="density"})}  
}

\details{A seats-votes curve (with various annotations) is produced
  with option \code{type="seatsVotes"}.

  A density plot of the vote shares is
  produced with \code{type="density"}. A bimodal density corresponds to
  an electoral system with a proliferation of safe seats for both
  parties, and a seats-votes curve that is relatively flat (or
  \dQuote{unresponsive}) in the neighborhood of average district-level
  vote shares of 50 percent.  The density is fitted using the defaults
  in the \code{\link{density}} function, but with the density constrained
  to fall to zero at the extremes of the data, via the \code{from} and
  \code{to} options to density.  A \code{\link{rug}} is added to the
  density plot.

  If \code{transform=TRUE}, the vote shares are transformed prior to
  plotting, so as to reduce the extent to which extreme vote shares
  close to zero or 1 determine the shape of the density (i.e., this
  option is available only for plots of \code{type="density"}).  The
  transformation is a sinusoidal function, a scaled
  \dQuote{log-odds/inverse-log-odds} function mapping from (0,1) to
  (0,1): i.e., \eqn{f(x) = g(k\cdot h(x))}{f(x) = g(kh(x))} where
  \eqn{h(x)} is the log-odds transformation \eqn{h(x) =
  \log(x/(1-x))}{h(x) = log[x/(1-x)]}, \eqn{k} is a scaling parameter
  set to \eqn{\sqrt{3}}{sqrt(3)}, and \eqn{g(x)} is the inverse-log-odds
  transformation \eqn{g(x) = \exp(x)/(1+\exp(x))}{g(x) =
  exp(x)/[1+exp(x)]}.  Note that this transformation is cosmetic, with
  the effect of assigning more of the graphing region to be devoted to
  marginal seats.
  
}

\value{
  The function performs the requested plots and exits silently with
  \code{invisible{NULL}}.
}

\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\seealso{\code{\link{density}}, \code{\link{rug}}}
\examples{
data(ca2006)
x <- ca2006$D/(ca2006$D+ca2006$R)
sv <- seatsVotes(x,
                 desc="Democratic Vote Shares, California 2006 congressional elections")

plot(sv)
plot(sv,type="density")
plot(sv,type="density",transform=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
