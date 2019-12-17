\name{readKH}
\alias{readKH}
\title{read roll call data in Poole-Rosenthal KH format}
\description{
  Creates a \code{rollcall} object from the flat file format for roll
  call data used by
  Keith Poole and Howard Rosenthal.
}

\usage{
readKH(file,
       dtl=NULL,
       yea=c(1,2,3),
       nay=c(4,5,6),
       missing=c(7,8,9),
       notInLegis=0,
       desc=NULL,
       debug=FALSE)
}

\arguments{
  \item{file}{string, name of a file or URL holding KH data}
  \item{dtl}{string, name of a file or URL holding KH \code{dtl} file
    (information about votes); default is \code{NULL}, indicating no
    \code{dtl} file} 
  \item{yea}{numeric, possibly a vector, code(s) for a Yea vote in the
    rollcall context (or a correct answer in the educational testing
    context). Default is \code{c(1,2,3)}, which corresponds to Yea,
    Paired Yea, and Announced Yea in Poole/Rosenthal data files.}
  \item{nay}{numeric, possibly a vector, code(s) for a Nay vote in the
    rollcall context (or an incorrect answer in the educational testing
    context).  Default is \code{c(4,5,6)}, which corresponds to
    Announced Nay, Paired Nay, and Nay in Poole/Rosenthal data files.}
  \item{missing}{numeric and/or \code{NA}, possible a vector, code(s)
    for missing data.  Default is \code{c(0,7,8,9,NA)}; the first four
    codes correspond to Not Yet a Member, Present (some Congresses),
    Present (some Congresses), and Not Voting.}
  \item{notInLegis}{numeric or \code{NA}, possibly a vector, code(s) for
    the legislator not being in the legislature when a particular roll
    call was recorded (e.g., deceased, retired, yet to be
    elected). Default is \code{0} for Poole/Rosenthal data files.}
  \item{desc}{string, describing the data, e.g., \code{82nd U.S. House
      of Representatives}; default is \code{NULL}}
  \item{debug}{logical, print debugging information for net connection}
}

\value{an object of class \code{\link{rollcall}}, with components created
  using the identifying information in the Poole/Rosenthal files.  If
  the function can not read the file (e.g., the user specified a URL and
  the machine is not connected to the Internet), the function fails with
  an error message (set \code{debug=TRUE} to help resolve these issues).
}

\details{Keith Poole and Howard Rosenthal have gathered an impressive
  collection of roll call data, spanning every roll call cast in
  the United States Congress.  This effort continues now as a real-time
  exercise, via a collaboration with Jeff Lewis (109th Congress
  onwards).  Nolan McCarty collaborated on the compilation of roll call
  data for the 102nd through 108th Congress.

  This function relies on some hard-coded features of Poole-Rosenthal
  flat files, and assumes that the \code{file} being supplied has the following
  structure (variable, start-end columns):
  \describe{
    \item{ICPSR legislator unique ID}{4-8}
    \item{ICPSR state ID}{9-10}
    \item{Congressional District}{11-12}
    \item{state name}{13-20}
    \item{party code}{21-23}
    \item{legislator name}{26-36}
    \item{roll-call voting record}{37 to end-of-record}
    }
  This function reads data files in that format, and creates a
  \code{\link{rollcall}}, for which there are useful methods such as
  \code{\link{summary.rollcall}}.  The \code{legis.data} component of the
  \code{\link{rollcall}} object is a \code{\link{data.frame}} which
  contains:
  \describe{
    \item{\code{state}}{a 2-character string abbreviation of each
      legislator' state}
    \item{\code{icpsrState}}{a 2-digit numeric code for each
      legislator's state, as used by the Inter-university Consortium for
      Political and Social Research (ICPSR)}
    \item{\code{cd}}{numeric, the number of each legislator's
      congressional district within each state; this is always 0 for
      members of the Senate}    
    \item{\code{icpsrLegis}}{a unique numeric identifier for each
      legislator assigned by the ICPSR, as corrected by Poole and
      Rosenthal.}
    \item{\code{partyName}}{character string, the name of each legislator's
      political party}
    \item{\code{party}}{numeric, code for each legislator's political
      party; see \url{http://legacy.voteview.com/PARTY3.HTM}}
  }
  The \code{\link{rownames}} attribute of this data frame is a
  concatenation of the legislators' names, party abbreviations (for
  Democrats and Republicans) and state, and (where appropriate), a
  district number; e.g., \code{Bonner (R AL-1)}.  This tag is also
  provided in the \code{legis.name} component of the returned rollcall
  object.

  Poole and Rosenthal also make \code{dtl} files available for
  Congresses 1 through 106. These files contain information about the
  votes themselves, in a multiple-line per vote \code{ascii} format, and
  reside in the \code{dtl} director of Poole's web site, e.g.,
  \url{https://legacy.voteview.com/k7ftp/dtl/102s.dtl} is the \code{dtl}
  file for the 102nd Senate.  The default is to presume that no such
  file exists.  When a \code{dtl} file is available, and is read, the
  \code{votes.data} attribute of the resulting \code{\link{rollcall}}
  object is a \code{\link{data.frame}} with one record per vote, with
  the following variables:
  \describe{
    \item{\code{date}}{vector of class \code{\link{Date}}, date of the
    rollcall, if available; otherwise \code{NULL}}
    \item{\code{description}}{vector of mode \code{character}, descriptive
      text}
  }
  The \code{dtl} files are presumed to have the date of the rollcall in
  the first line of text for each roll call, and lines 3 onwards contain
  descriptive text.
  
  Finally, note also that the Poole/Rosenthal data sets often include
  the U.S. President as a pseudo-legislator, adding the announced
  positions of a president or the administration to the roll call
  matrix. This adds an extra \dQuote{legislator} to the data set and can
  sometimes produce surprising results (e.g., a U.S. Senate of 101
  senators), and a \dQuote{legislator} with a surprisingly low party loyalty
  score (since the President/administration only announces positions on
  a relatively small fraction of all Congressional roll calls).   
}

\references{Poole, Keith and Howard Rosenthal.  1997.  \emph{Congress: A
  Political-Economic History of Roll Call Voting}.  New York: Oxford
  University Press.  

Poole, Keith.  \url{http://legacy.voteview.com}

Rosenthal, Howard L. and Keith T. Poole.  \emph{United States Congressional
Roll Call Voting Records, 1789-1990: Reformatted Data [computer
file].} 2nd ICPSR release. Pittsburgh, PA: Howard L. Rosenthal and Keith
T. Poole, Carnegie Mellon University, Graduate School of Industrial
Administration [producers], 1991. Ann Arbor, MI: Inter-university
Consortium for Political and Social Research [distributor], 2000.  
\url{http://www.icpsr.umich.edu/icpsrweb/ICPSR/studies/09822}
}

\seealso{\code{\link{rollcall}}
}

\examples{
\dontrun{
h107 <- readKH("https://voteview.com/static/data/out/votes/H107_votes.ord",
                desc="107th U.S. House of Representatives")

s107 <- readKH("https://voteview.com/static/data/out/votes/S107_votes.ord",
                desc="107th U.S. Senate")
}
}
\author{Simon Jackman \email{simon.jackman@sydney.edu.au}}
\keyword{datasets}

