\name{maxdegree-term}
\alias{maxdegree}
\docType{methods}
\title{Maximum Degree Term for the 'ergm' Package}
\description{
  This term counts the number of nodes in the network with a given degree or less.
}
%\usage{
%terms.ergm()
\section{The Maximum Degree Term}{
 \describe{
  \item{\code{maxdegree(maxdeg, by=NULL)}}{\emph{Minimum Degree:}
    The \code{maxdeg} argument is an integer. This
term adds one
    network statistic to the model being the 
number of nodes in the network of at most degree \code{maxed}.
That is, the statistic equals the number of nodes in the
network with \code{maxdeg} or less edges.
The optional argument
    \code{by} is a character string giving the name of an
attribute in the
    network's vertex attribute list.
    If this is specified, 
    then degrees are calculated using the subnetwork consisting
of only
    edges whose endpoints have the same value of the \code{by}
attribute.
    This term can only be used with undirected networks.
   }
  }
}
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

\details{
     
 This is a trivial modification of the \code{mindegree} term in the \code{ergm.userterms} package and the ergm.useterms tutorials. 
}
%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\author{
David Hunter, Steve Goodreau, Mark Handcock, Carter Butts 
}

\references{
Hunter DR, Goodreau SM, Handcock MS (2013).
{\pkg{ergm.userterms}: {A} Template Package for 
Extending \pkg{statnet}}, Journal of Statistical
Software 52(2), 1-25,
URL http://www.jstatsoft.org/v52/i02/.
}
\seealso{statnet, network, ergm, ergm-terms}
\examples{
data(florentine)
summary(flomarriage~maxdegree(3))
summary(flomarriage~maxdegree(1,by="priorates"))
\donttest{
fit <- ergm(flomarriage~edges+maxdegree(1,by="priorates"))
summary(fit)
}
}
\keyword{models}
