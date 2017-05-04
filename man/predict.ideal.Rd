\name{predict.ideal}
\alias{predict.ideal}
\alias{print.predict.ideal}
\title{predicted probabilities from an ideal object}

\description{Compute predicted probabilities from an \code{\link{ideal}}
  object.  This predict method uses the posterior mean values of \eqn{x}
  and \eqn{\beta}{beta} to make predictions.}

\usage{
\method{predict}{ideal}(object,
                        cutoff=.5,
                        burnin=NULL,
                        ...)

\method{print}{predict.ideal}(x,digits=2,...)

}
\arguments{
  \item{object}{an object of class \code{ideal} (produced by
    \code{\link{ideal}}) with item parameters (beta) stored; i.e.,
    \code{store.item=TRUE} was set when the \code{ideal} object was
    fitted} 
  \item{cutoff}{numeric, a value between 0 and 1, the threshold to be
    used for classifying predicted probabilities of a Yea 
    votes as predicted Yea and Nay votes.}
  \item{burnin}{of the recorded MCMC samples, how many to discard as
    burnin?  Default is \code{NULL}, in which case the value of
    \code{burnin} in the \code{\link{ideal}} object is used.}
  \item{x}{object of class \code{predict.ideal}}
  \item{digits}{number of digits in printed object}
  \item{...}{further arguments passed to or from other methods.}  
}

\details{
  Predicted probabilities are computed using the mean of the posterior
  density of
  of \eqn{x} (ideal points, or latent ability) and \eqn{\beta} (bill or
  item parameters).  The percentage correctly predicted
  are determined by counting the percentages of votes with predicted
  probabilities of a Yea vote greater than or equal to the \code{cutoff} as the
  threshold.
}

\value{
  An object of class \code{predict.ideal}, containing:
  \item{pred.probs}{the calculated predicted probability for each
    legislator for each vote.}
  \item{prediction}{the calculated prediction (0 or 1) for each
    legislator for each vote.}
  \item{correct}{for each legislator for each vote, whether the
    prediction was correct.}
  \item{legis.percent}{for each legislator, the percent of votes
    correctly predicted.}
  \item{vote.percent}{for each vote, the percent correctly predicted.}
  \item{yea.percent}{the percent of yea votes correctly predicted.}
  \item{nay.percent}{the percent of nay votes correctly predicted.}
  \item{party.percent}{the average value of the percent correctly
    predicted by legislator, separated by party, if party information
    exists in the \code{rollcall} object used for \code{ideal}.  If no
    party information is available, \code{party.percent = NULL}.}
  \item{overall.percent}{the total percent of votes correctly
    predicted.}
  \item{ideal}{the name of the \code{\link{ideal}} object, which can be
    later \code{\link{eval}}uated}
  \item{desc}{string, the descriptive text from the
  \code{\link{rollcall}} object passed to \code{\link{ideal}}}
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

\seealso{\code{\link{ideal}}, \code{\link{summary.ideal}}, \code{\link{plot.predict.ideal}}}

\examples{
data(s109)

f <- system.file("extdata","id1.rda",package="pscl")
load(f)
phat <- predict(id1)
phat         ## print method
}

\keyword{classes}
