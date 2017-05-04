\name{simpi}
\alias{simpi}

\title{Monte Carlo estimate of pi (3.14159265...)}
\description{
  Monte Carlo estimation of pi
}
\usage{
simpi(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{integer, number of Monte Carlo samples, defaults to 1000}
}

\details{
  A crude Monte Carlo estimate of \eqn{\pi}{pi} can be formed as
  follows.  Sample from the unit square many times (i.e., each sample is
  formed with two independent draws from a uniform density on the unit
  interval).  Compute the proportion \eqn{p}{p} of sampled points that
  lie inside a unit circle centered on the origin; such points
  \eqn{(x,y)}{(x,y)} have distance from the origin \eqn{d = \sqrt{x^2 +
  y^2}}{d=sqrt(x^2 + y^2)} less than 1.  Four times \eqn{p}{p} is a
  Monte Carlo estimate of \eqn{\pi}{pi}.  This function is a wrapper to
  a simple C function, bringing noticable speed gains and memory
  efficiencies over implementations in native R.

  Contrast this Monte Carlo method with Buffon's needle and refinements thereof (see
  the discussion in Ripley (1987, 193ff). 
}
  
\value{the Monte Carlo estimate of \eqn{\pi}{pi}}

\references{Ripley, Brain D. 1987 [2006].  \emph{Stochastic Simulation}.
  Wiley: Hoboken, New Jersey.}

\author{Simon Jackman \email{jackman@stanford.edu}}

\examples{
seed <- round(pi*10000)  ## hah hah hah
m <- 6
z <- rep(NA,m)
lim <- rep(NA,m)
for(i in 1:m){
  cat(paste("simulation for ",i,"\n"))
  set.seed(seed)
  timings <- system.time(z[i] <- simpi(10^i))
  print(timings)
  cat("\n")
  lim[i] <- qbinom(prob=pi/4,size=10^i,.975)/10^i * 4
}

## convert to squared error
z <-(z - pi)^2
lim <- (lim - pi)^2

plot(x=1:m,
     y=z,
     type="b",
     pch=16,
     log="y",
     axes=FALSE,
     ylim=range(z,lim),
     xlab="Monte Carlo Samples",
     ylab="Log Squared Error")
lines(1:m,lim,col="blue",type="b",pch=1)
legend(x="topright",
       legend=c("95\% bound",
         "Realized"),
       pch=c(1,16),
       lty=c(1,1),
       col=c("blue","black"),
       bty="n")
axis(1,at=1:m,
     labels=c(expression(10^{1}),
       expression(10^{2}),
       expression(10^{3}),
       expression(10^{4}),
       expression(10^{5}),
       expression(10^{6})))
axis(2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
