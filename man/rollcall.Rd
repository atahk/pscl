\name{rollcall}
\alias{rollcall}
\title{create an object of class rollcall}
\description{
  Create a \code{rollcall} object, used for the analysis of legislative
  voting or, equivalently, item-response modeling of binary data
  produced by standardized tests, etc.
}

\usage{
rollcall(data,
         yea=1, nay=0, missing=NA, notInLegis=9,
         legis.names=NULL, vote.names=NULL,
         legis.data=NULL, vote.data=NULL,
         desc=NULL, source=NULL)
}

\arguments{

  \item{data}{voting decisions (for roll calls), or test results (for
    IRT).  Can be in one of two forms.  First, \code{data} may be a
    \code{\link{matrix}}, with rows corresponding to legislators
    (subjects) and columns to roll calls (test items).  \code{data} can
    also be a \code{\link{list}} with an element named \code{votes}
    containing the matrix described above.}
  
  \item{yea}{numeric, possibly a vector, code(s) for a Yea vote in the
    rollcall context, or a correct answer in the educational testing
    context. Default is 1.}
  
  \item{nay}{numeric, possibly a vector, code(s) for a Nay vote in the
    rollcall context, or an incorrect answer in the educational testing
    context.  Default is 0.}
  
  \item{missing}{numeric or \code{NA}, possibly a vector, code(s) for
    missing data.  Default is \code{NA}.}
  
  \item{notInLegis}{numeric or \code{NA}, possibly a vector, code(s) for
    the legislator not being in the legislature when a particular roll
    call was recorded (e.g., deceased, retired, yet to be elected).}
  
  \item{legis.names}{a vector of names of the legislators or
    individuals.  If \code{data} is a \code{list} or \code{data.frame}
    and has a component named \code{legis.names}, then this will be
    used.  Names will be generated if not supplied, or if there are
    fewer unique names supplied than legislators/subjects (rows of the
    roll call matrix).}
  
  \item{vote.names}{a vector of names or labels for the votes or items.
    If \code{data} is a \code{list} or \code{data.frame} and has a
    component named \code{vote.names}, then this will be used.  Names
    will be generated if not supplied, or if there are fewer unique
    names supplied than votes/test-items (columns of the roll call
    matrix).}
  
  \item{legis.data}{a \code{\link{matrix}} or \code{\link{data.frame}}
    containing covariates specific to each legislator/test-taker; e.g.,
    party affiliation, district-level covariates.  If this object does
    not have the same number of rows as \code{data}, an error is
    returned.}
  
  \item{vote.data}{a \code{\link{matrix}} or \code{\link{data.frame}}
    containing covariates specific to each roll call vote or test item;
    e.g., a timestamp, the bill sponsor, descriptive text indicating
    the type of vote.  If this object does not have the same number of
    row as the number of columns in \code{data}, an error is returned.}
  
  \item{desc}{character, a string providing an (optional) description of
    the data being used.  If \code{data} is a list and contains an
    element named \code{desc}, then this will be used.}
  
  \item{source}{character, a string providing an (optional) description
    of where the roll call data originated (e.g., a URL or a short-form
    reference).  Used in print and summary methods.} 
}

\details{See below for methods that operate on objects of class
  \code{rollcall}.} 

\value{
  An object of class \code{rollcall}, a list with the following components:
  
  \item{votes}{a \code{\link{matrix}} containing voting decisions, with
    rows corresponding to legislators (test subjects) and columns to
    roll call votes (test items). Legislators (test subjects)
    and items (or votes) have been labeled in the
    \code{\link{dimnames}} attribute of this matrix, using
    the \code{legis.names} and/or \code{vote.names} arguments, respectively.}

  \item{codes}{a \code{\link{list}} with named components \code{yea},
    \code{nay}, \code{notInLegis} and \code{missing}, each component a
    numeric vector (possibly of length 1 and possibly \code{NA}),
    indicating how entries in the \code{votes} component of the
    \code{rollcall} object should be considered.  This list simply gathers
    up the values in the \code{yea}, \code{nay}, \code{notInLegis} and
    \code{missing} arguments passed to the function.}
  
  \item{n}{numeric, number of legislators, equal to \code{dim(votes)[1]}}

  \item{m}{numeric, number of votes, equal to \code{dim(votes)[2]}}

  \item{legis.data}{user-supplied data on legislators/test-subjects.}
  
  \item{vote.data}{user-supplied data on rollcall votes/test-items.}
  
  \item{desc}{any user-supplied description.  If no description was provided,
    defaults \code{desc} defaults to \code{NULL}.}
  
  \item{source}{any user-supplied source information (e.g., a url or a
    short-form reference).  If no description is provided, \code{source}
    defaults to \code{NULL}.}
  
}

\seealso{
  \code{\link{readKH}} for creating objects from files (possibly over
  the web), in the format used for data from the United States Congress
  used by Keith Poole and Howard Rosenthal (and others).

  \code{\link{summary.rollcall}}, \code{\link{ideal}} for model fitting.
}

\examples{
## generate some fake roll call data
set.seed(314159265)
fakeData <- matrix(sample(x=c(0,1),size=5000,replace=TRUE),
                   50,100)
rc <- rollcall(fakeData)
is(rc,"rollcall")        ## TRUE
rc                       ## print the rollcall object on screen

data(sc9497)             ## Supreme Court example data
rc <- rollcall(data=sc9497$votes,
               legis.names=sc9497$legis.names,
               desc=sc9497$desc)
summary(rc,verbose=TRUE)                      

\dontrun{
## s107
## could use readKH for this
dat <- readLines("sen107kh.ord")
dat <- substring(dat,37)
mat <- matrix(NA,ncol=nchar(dat[1]),nrow=length(dat))
for(i in 1:103){
  mat[i,] <- as.numeric(unlist(strsplit(dat[i],
                                        split=character(0))))
}

s107 <- rollcall(mat,
                 yea=c(1,2,3),
                 nay=c(4,5,6),
                 missing=c(7,8,9),
                 notInLegis=0,
                 desc="107th U.S. Senate",
                 source="http://voteview.ucsd.edu")
summary(s107)
}
}
\keyword{manip}



