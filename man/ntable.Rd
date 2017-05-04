\name{ntable}
\alias{ntable}
\title{nicely formatted tables}
\description{
  Nicely formatted tables, with row or column marginals etc.
}
\usage{
ntable(x,y=NULL,
       percent=1,digits=2,
       row=FALSE,col=FALSE)
}

\arguments{
  \item{x}{vector or \code{\link{factor}}}
  \item{y}{vector of \code{\link{factor}}}
  \item{percent}{integer, 1 for row percentages (default), 2 for column
    percentages}
  \item{digits}{integer, digits to print after decimal place (default is 2)}
  \item{row}{logical, if \code{TRUE}, print row marginals}
  \item{col}{logical, if \code{TRUE}, print column marginals}
}

\details{A wrapper function to \code{\link{prop.table}} that produces prettier
  looking results.
}

\value{
  nothing returned; the function prints the table and
  exits silently.
}
\seealso{\code{\link{prop.table}}, \code{\link{table}}}

\author{Jim Fearon \email{jfearon@stanford.edu}}

\examples{
data(bioChemists)
attach(bioChemists)
ntable(fem)
ntable(fem,mar,row=TRUE)
ntable(fem,mar,per=2,col=TRUE)
ntable(fem,mar,per=2,row=TRUE,col=TRUE)
}
\keyword{print}
