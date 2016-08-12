\name{gwbnsp-term}
\alias{gwbnsp}
\alias{gwb1nsp}
\alias{gwb2nsp}
\docType{methods}
\title{Geometrically weighted bipartite non-edgewise shared partner terms (gwb1nsp and gwb2nsp)}
\description{
  Geometrically weighted non-edgewise shared partner distribution for first and second modes of bipartite networks. These terms are closely related to the \code{gwnsp} term in the \pkg{ergm} package, with adaptations to work for bipartite networks and distinguish between nsp's on each mode
}
%\usage{
%terms.ergm()

%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

\details{

\code{gwb1nsp(alpha=0, fixed=TRUE) (binary) (undirected) (curved)}

Geometrically weighted nonedgewise shared partner distribution for the first mode in a bipartite (aka two-mode) network: This term is just like gwnsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with weight parameter alpha for the first mode in a bipartite network. The argument fixed indicates whether the scale parameter lambda is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006). Currently, the fixed argument must be set to TRUE, which means the weight parameter alpha is not separately estimated and thus the model is not a CEF model. This term can only be used with undirected bipartite networks.

\code{gwb2nsp(alpha=0, fixed=TRUE) (binary) (undirected) (curved)}

Geometrically weighted nonedgewise shared partner distribution for the second mode in a bipartite (aka two-mode) network: This term is just like gwnsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with weight parameter alpha for the second mode in a bipartite network. The argument fixed indicates whether the scale parameter lambda is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006). Currently, the fixed argument must be set to TRUE, which means the weight parameter alpha is not separately estimated and thus the model is not a CEF model. This term can only be used with undirected bipartite networks.

Note that allthough the default is \code{fixed=FALSE} the terms are not yet able to handle a non-fixed decay term so it must be set to TRUE.
}

\seealso{
  See also \code{\link[ergm]{gwnsp}} term. 
}


%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\references{
Hunter, D. R. and M. S. Handcock (2006). Inference in curved exponential family models for networks. Journal of Computational and Graphical Statistics, 15: 565-583.

Wang, P., Sharpe, K., Robins, G. L., and Pattison, P. E. (2009). Exponential random graph (p*) models for affiliation networks. Social Networks. 31(1): 12-25.

}

\author{
Matthew Hamilton  \email{matthewlhamilton@gmail.com}
}

\examples{
require(ergm.terms.contrib)
# create example data
bel <- data.frame(actor = c("a","a","a","b","b","b","b","c","c","c","c"), 
                  event = c("e1","e2","e3","e3","e1","e4","e5","e3","e1","e4","e6"))
bnet <- network(bel,bipartite=3,directed=F)

# how many actor (first mode) dyads have ANY shared partners?
summary(bnet ~ gwb1nsp(0,TRUE))

# calculate for 2nd mode, using non-extreme alphas 
# to weight each additional shared partner less than the one before
summary(bnet ~ gwb2nsp(.5,TRUE))

}
\keyword{models}
