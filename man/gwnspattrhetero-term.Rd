\name{gwnspattrhetero-term}
\alias{gwnspattrhetero}
\docType{methods}
\title{Geometrically weighted nonedgewise shared partner term for dyads with heterogeneous attributes}
\description{
  Geometrically weighted nonedgewise shared partner distribution for dyads with heterogeneous attributes. This term is closely related to the \code{gwnsp} term in the \pkg{ergm} package.
}
%\usage{
%terms.ergm()
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}
\details{
\code{gwnspattrhetero(attrname=NULL,pexcl=NULL,alpha=NULL) (curved)}
Geometrically weighted nonedgewise shared partner distribution for dyads with heterogeneous attributes: Similar to the gwnsp term, this term adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with weight parameter alpha; however the statistic only counts nonedgewise shared partners of dyads with different values of the attribute specified by the attrname argument. The argument pexcl allows one attribute value to to excluded. Currently, the pexcl argument is required; if no values need to be excluded, this argument can be set to 99 or some number higher than the number of values of the attribute. This term can only be used with directed networks.
}
\seealso{
  See also \code{\link[ergm]{gwnsp}} term. 
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
net <- network.initialize(5,directed=T)
net%v%"pendants" <- c("zgreen","red","black", "zgreen", "red")
net[1,2] <- net[1,3] <- net[3,4] <- net[4,2] <- net[3,5] <- net[2,5] <- 1

plot(net, vertex.col = as.factor(net%v%"pendants"), displaylabels = T)

# how many non-connected dyads with different attribute values have ANY shared partners?
summary(net ~ dgwnspattr3("pendants",99,0))
# 1->X->5; 4->X->5; 3->X->2

# how many shared partners are there for non-connected dyads with different attribute values?
summary(net ~ dgwnspattr3("pendants",99,100))
# 1->X->5 (2 shared partners); 4->X->5 (1 shared partner); 3->X->2 (1 shared partner)

# non-extreme alphas weight each additional shared partner less than one before
summary(net ~ dgwnspattr3("pendants",99,1))
}

\keyword{models}
