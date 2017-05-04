\name{summary.rollcall}
\alias{summary.rollcall}
\alias{print.summary.rollcall}
\title{summarize a rollcall object}

\description{
  Provides a summary of the information about votes, legislators, etc in
  a \code{\link{rollcall}} object. 
}

\usage{
\method{summary}{rollcall}(object,
                 dropList=NULL,  
                 verbose=FALSE,debug=FALSE,...)

\method{print}{summary.rollcall}(x, digits=1, ...)
}

\arguments{
  \item{object}{an \code{\link{rollcall}} object.}
  \item{dropList}{a \code{\link{list}} or \code{\link{alist}}, listing voting
    decisions, legislators and/or votes to be dropped 
    from the summary; see \code{\link{dropRollCall}} for details.}
  \item{verbose}{logical, if \code{TRUE}, compute legislator-specific
    and vote-specific Yea/Nay/\code{NA} summaries}
  \item{debug}{logical, if \code{TRUE}, print messages to console during
  processing of the \code{rollcall} object}
  \item{x}{an object of class \code{summary.rollcall}}
  \item{digits}{number of decimal places in printed display}
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  An object of class \code{summary.rollcall} with the following elements
  (depending on the logical flag \code{verbose}):

  \item{n}{number of legislators in the \code{\link{rollcall}} object,
    after processing the \code{dropList}}

  \item{m}{number of roll call votes in the \code{\link{rollcall}}
    object, after processing the \code{dropList}}
  
  \item{codes}{a \code{\link{list}} that describes how the voting
    decisions in the \code{\link{rollcall}} matrix (\code{object$votes})
    map into \dQuote{Yea} and \dQuote{Nay} etc, after processing the
    \code{dropList}; see \code{\link{rollcall}} for more details}
  
  \item{allVotes}{a matrix containing a tabular breakdown of all votes
    in the \code{\link{rollcall}} matrix (\code{object$votes}), after
    processing the \code{dropList}}

  \item{partyTab}{a tabular breakdown of the legislators' party
    affiliations, after processing the \code{dropList}, and only if
    party affiliations are supplied as\code{object$legis.data$party};
    see \code{\link{rollcall}} for details}
  
  \item{lopSided}{a tabular summary of the frequency of lop-sided roll
      call votes in the \code{\link{rollcall}} object, again, after
      processing the \code{dropList}}
    
    \item{legisTab}{a tabular summary of each legislators' voting
      history}

    \item{partyLoyalty}{the proportion of times that each legislator
    votes the way that a majority of his or her fellow partisans did,
    provided party affiliations are available}

  \item{voteTab}{a tabular summary of each rollcall's votes}

  \item{call}{the \code{\link[=match.call]{matched call}} used to invoke
    \code{summary.rollcall}}
    
}

\seealso{\code{\link{rollcall}}}

\examples{
set.seed(314159265)
fakeData <- matrix(sample(x=c(0,1),size=1000,replace=TRUE),
                   10,100)
rc <- rollcall(fakeData)
rc

data(sc9497)
rc <- rollcall(sc9497)
summary(rc)

data(s109)
summary(s109)
summary(s109,verbose=TRUE)
}

\keyword{classes}

