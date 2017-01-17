\name{gwdspattr-term}
\alias{gwdspattr}
\docType{methods}
\title{Geometrically weighted dyadwise shared partner terms with dyad attributes}
\description{
  Geometrically weighted dyadwise shared partner distribution for which attibutes of dyads are specified. This term is closely related to the \code{gwdsp} term in the \pkg{ergm} package.
}
%\usage{
%terms.ergm()
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}
\details{
\code{gwdspattr(attrname=NULL,p1keep=NULL,p2keep=NULL,alpha=NULL) (curved)}
Geometrically weighted nonedgewise shared partner distribution for the first mode in a bipartite (aka two-mode) network: This term is just like gwnsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with weight parameter alpha for the first mode in a bipartite network. The argument fixed indicates whether the scale parameter lambda is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006). Currently, the fixed argument must be set to TRUE, which means the weight parameter alpha is not separately estimated and thus the model is not a CEF model. This term can only be used with undirected bipartite networks.
}
\seealso{
  See also \code{\link[ergm]{gwdsp}} term. 
}
%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\references{
Hunter, D. R. and M. S. Handcock (2006). Inference in curved exponential family models for networks. Journal of Computational and Graphical Statistics, 15: 565-583.
}
\author{
Matthew Hamilton  \email{matthewlhamilton@gmail.com}
}
\examples{
require(ergm.userterms.contrib)
require(sna)

# create example data
net <- network.initialize(5,directed=F)
net%v%"pendants" <- c("zgreen","red","black", "zgreen", "red")
net[1,2] <- net[1,3] <- net[3,4] <- net[4,2] <- net[3,5] <- net[2,5] <- 1

# how many black-red dyads have ANY shared partners?
summary(net ~ gwdspattr("pendants",1,2,0))

# how many shared partners are there for all black-red dyads?
summary(net ~ gwdspattr("pendants",1,2,100))

# non-extreme alphas weight each additional shared partner less than one before
summary(net ~ gwdspattr("pendants",1,2,1))
}

\keyword{models}
