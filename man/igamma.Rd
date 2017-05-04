\name{igamma}
\alias{igamma}
\alias{densigamma}
\alias{pigamma}
\alias{qigamma}
\alias{rigamma}
\alias{igammaHDR}
\title{inverse-Gamma distribution}

\description{Density, distribution function, quantile function, and highest density region calculation for the inverse-Gamma distribution with parameters \code{alpha} and \code{beta}.}

\usage{
	densigamma(x,alpha,beta)
  	pigamma(q,alpha,beta)
	qigamma(p,alpha,beta)
	rigamma(n,alpha,beta)
	igammaHDR(alpha,beta,content=.95,debug=FALSE)
}

\arguments{
	\item{x,q}{vector of quantiles}
	\item{p}{vector of probabilities}
	\item{n}{number of random samples in \code{rigamma}} 
	\item{alpha,beta}{rate and shape parameters of the
	inverse-Gamma density, both positive}
	\item{content}{scalar, 0 < \code{content} < 1, volume of
	highest density region}
	\item{debug}{logical; if TRUE, debugging information from the
	search for the HDR is printed} 
    }
    
\details{
  The inverse-Gamma density arises frequently in Bayesian
  analysis of normal data, as the (marginal) conjugate prior for
  the unknown variance parameter.  The inverse-Gamma density
  for \eqn{x>0} with parameters \eqn{\alpha>0} and \eqn{\beta>0}
  is
  \deqn{f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{-\alpha-1}
    \exp(-\beta/x)}{%
    (beta^alpha)/Gamma(alpha) x^(-alpha-1) exp(-beta/x)
  }
  where \eqn{\Gamma(x)} is the \code{\link{gamma}} function
  \deqn{\Gamma(a) = \int_0^\infty t^{a-1} \exp(-t) dt}{%
    Gamma(a) = int_0^infty t^(a-1) exp(-t) dt
  }
  and so ensures \eqn{f(x)} integrates to one.  The
  inverse-Gamma density has a mean at \eqn{\beta/(\alpha-1)}{beta/(alpha-1)} for
  \eqn{\alpha>1}{alpha>1} and has variance \eqn{\beta^2/((\alpha-1)^2
    (\alpha-2))}{beta^2/((alpha-1)^2 (alpha-2))} for
  \eqn{\alpha>2}{alpha>2}.  The inverse-Gamma density has a unique mode at
  \eqn{\beta/(\alpha+1)}{beta/(alpha+1)}.  
  
  The evaluation of the density,
  cumulative distribution function and quantiles is done by
  calls to the \code{dgamma}, \code{pgamma} and \code{igamma}
  functions, with the arguments appropriately transformed.  That
  is, note
  that if \eqn{x \sim IG(\alpha,\beta)}{x ~ IG(alpha,beta} then \eqn{1/x \sim
    G(\alpha,\beta)}{1/x ~ G(alpha,beta)}.  
  
  \emph{Highest Density Regions}. In general, suppose \eqn{x}
  has a density \eqn{f(x)}, where \eqn{x \in \Theta}.  Then a
  highest density region (HDR) for \eqn{x} with content \eqn{p
    \in (0,1]} is a region (or set of regions) \eqn{\mathcal{Q}
    \subseteq \Theta} such that:
  \deqn{\int_\mathcal{Q} f(x) dx = p}{%
    int_Q f(x) dx = p
  }
  and
  \deqn{f(x) > f(x^*) \, \forall\ x \in \mathcal{Q},
    x^* \not\in \mathcal{Q}.}{%
    f(x) > f(x*) for all x in Q and all x* not in Q.
  }
  For a continuous, unimodal
  density defined with respect to a single parameter (like the
  inverse-Gamma case considered here with parameters \eqn{0 <
    \alpha < \infty, \,\, 0 < \beta < \infty}), a HDR region \eqn{Q}
  of content \eqn{p} (with \eqn{0 < p < 1}) is a unique, closed
  interval on the real half-line.
  
  This function uses numerical methods to solve for the
  boundaries of a HDR with \code{content} \eqn{p} for the
  inverse-Gamma density, via repeated calls the functions
  \code{densigamma}, \code{pigamma} and
  \code{qigamma}. In particular, the function \code{\link{uniroot}} is used to
  find points \eqn{v} and \eqn{w} such that \deqn{f(v) = f(w)} subject
  to the constraint
  \deqn{\int_v^w f(x; \alpha, \beta) d\theta = p.}{%
    int_v^w f(x; alpha, beta) d theta = p.
    }
}

\value{\code{densigamma} gives the density, \code{pigamma} the
  distribution function, \code{qigamma} the quantile function,
  \code{rigamma} generates random samples, and \code{igammaHDR} gives
  the lower and upper limits of the HDR, as defined above (\code{NA}s if
  the optimization is not successful).
}

\note{The \code{densigamma} is named so as not to conflict with the
  \code{\link{digamma}} function in the R \code{\link{base}} package
  (the derivative of the \code{gamma} function).
}

\author{Simon Jackman \email{jackman@stanford.edu}}

\seealso{\code{\link{gamma}}, \code{\link{dgamma}},
\code{\link{pgamma}}, \code{\link{qgamma}}, \code{\link{uniroot}}}

\examples{
alpha <- 4
beta <- 30
summary(rigamma(n=1000,alpha,beta))

xseq <- seq(.1,30,by=.1)
fx <- densigamma(xseq,alpha,beta)
plot(xseq,fx,type="n",
     xlab="x",
     ylab="f(x)",
     ylim=c(0,1.01*max(fx)),
     yaxs="i",
     axes=FALSE)
axis(1)
title(substitute(list(alpha==a,beta==b),list(a=alpha,b=beta)))
q <- igammaHDR(alpha,beta,debug=TRUE)
xlo <- which.min(abs(q[1]-xseq))
xup <- which.min(abs(q[2]-xseq))
plotZero <- par()$usr[3]
polygon(x=xseq[c(xlo,xlo:xup,xup:xlo)],
        y=c(plotZero,
          fx[xlo:xup],
          rep(plotZero,length(xlo:xup))),
        border=FALSE,
        col=gray(.45))
lines(xseq,fx,lwd=1.25)


\dontrun{
alpha <- beta <- .1
xseq <- exp(seq(-7,30,length=1001))
fx <- densigamma(xseq,alpha,beta)
plot(xseq,fx,
     log="xy",
     type="l",
     ylim=c(min(fx),1.01*max(fx)),
     yaxs="i",
     xlab="x, log scale",
     ylab="f(x), log scale",
     axes=FALSE)
axis(1)

title(substitute(list(alpha==a,beta==b),list(a=alpha,b=beta)))
q <- igammaHDR(alpha,beta,debug=TRUE)
xlo <- which.min(abs(q[1]-xseq))
xup <- which.min(abs(q[2]-xseq))
plotZero <- min(fx)
polygon(x=xseq[c(xlo,xlo:xup,xup:xlo)],
        y=c(plotZero,
          fx[xlo:xup],
          rep(plotZero,length(xlo:xup))),
        border=FALSE,
        col=gray(.45))
lines(xseq,fx,lwd=1.25)
}
}
\keyword{distribution}
