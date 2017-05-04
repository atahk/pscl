\name{prussian}
\alias{prussian}
\docType{data}

\title{Prussian army horse kick data}

\description{
  Deaths by year, by corp, from horse kicks.
}

\usage{data(prussian)}

\format{
  A data frame with 280 observations on the following 3 variables.
  \describe{
    \item{\code{y}}{a numeric vector, count of deaths}
    \item{\code{year}}{a numeric vector, 18XX, year of observation}
    \item{\code{corp}}{a \code{\link{factor}}, corp of Prussian Army generating
      observation}
  }
}

\source{
  von Bortkiewicz, L. 1898. \emph{Das Gesetz der Kleinen Zahlen.} Leipzig: Teubner. 
}

\examples{
data(prussian)
corpP <- glm(y ~ corp, family=poisson,data=prussian)
summary(corpP)
}
\keyword{datasets}
