\name{esp-term}
\alias{desp}
\alias{ddsp}
\alias{dnsp}
\alias{InitErgmTerm.desp}
\alias{InitErgmTerm.ddsp}
\alias{InitErgmTerm.dnsp}
\alias{InitErgmTerm.qwesp}

\alias{qwesp}
\alias{ergm-terms}

\alias{InitErgmTerm.dgwesp}
\alias{InitErgmTerm.dgwdsp}
\alias{InitErgmTerm.dgwnsp}
\alias{dgwesp}
\alias{dgwdsp}
\alias{dgwnsp}
\docType{methods}
\title{Extended Shared Partner Terms for directed networks for the 'ergm' Package}
\description{
  These terms extend (and replace) the existing terms within the \code{\link{ergm}} package dealing with edgewise shared partners; they differ from the conventional terms primarily in supporting multiple types of shared partner definitions in the directed case. 
}
%\usage{
%terms.ergm()
\section{ERGM Terms}{
 \describe{
  \item{\code{desp(d, type="OTP")}}{ (binary) \emph{Edgewise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{edges} (rather than dyads) in the network with exactly \code{d[i]} shared partners. This term can be used with directed and undirected networks.  In the directed case, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to be counted (see above for type codes).  By default, outgoing two-paths are employed.
   }
   
  \item{\code{dgwesp(alpha=0, fixed=FALSE, cutoff=30, type="OTP")}}{ (binary) \emph{Geometrically weighted edgewise shared partner distribution:} This term adds a statistic equal to the geometrically weighted \emph{edgewise} (not dyadwise) shared partner distribution with weight parameter \code{alpha}. The optional argument \code{fixed} indicates whether the scale parameter \code{lambda} is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006).  The default is \code{FALSE}, which means the scale parameter is not fixed and thus the model is a CEF model. This term can be used with directed and undirected networks. For directed networks, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see above for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying ESP terms to use in computing the statistics to reduce the computational burden.
   }
  
   \item{\code{dgwnsp(alpha=0, fixed=FALSE, cutoff=30, type="OTP")}}{ (binary) \emph{Geometrically weighted non-edgewise shared partner distribution:} This term is just like gwesp and gwdsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with weight parameter alpha. The optional argument \code{fixed} indicates whether the scale parameter \code{lambda} is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006).  The default is \code{FALSE}, which means the scale parameter is not fixed and thus the model is a CEF model. This term can be used with directed and undirected networks. For directed networks, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see above for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying ESP terms to use in computing the statistics to reduce the computational burden.
   }
   
   \item{\code{dgwdsp(alpha=0, fixed=FALSE, cutoff=30, type="OTP")}}{ (binary) \emph{Geometrically weighted dyadwise shared partner distribution:} This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared partner distribution with weight parameter alpha > 0. Note that the GWDSP statistic is equal to the sum of GWNSP plus GWESP. The optional argument \code{fixed} indicates whether the scale parameter \code{lambda} is to be fit as a curved exponential-family model (see Hunter and Handcock, 2006).  The default is \code{FALSE}, which means the scale parameter is not fixed and thus the model is a CEF model. This term can be used with directed and undirected networks. For directed networks, multiple shared partner definitions are possible; the \code{type} argument may be used to select the type of shared partner to employ (see above for definitions).  By default, outgoing two-paths are employed.  The optional argument \code{cutoff} sets the number of underlying ESP terms to use in computing the statistics to reduce the computational burden.
   }
   
  }
}
%}
%\arguments{
%\item{\dots}{Specific to the model term.}
%}

\details{
While there is only one shared partner configuration in the undirected case, nine distinct configurations are possible for directed graphs.  Currently, edgewise shared partner terms may be defined with respect to five of these configurations; they are defined here as follows (using terminology from Butts (2008) and the \code{relevent} package):

  \describe{
    \item{Outgoing Two-path (OTP)}{vertex \eqn{k} is an OTP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k \to j}{i->k->j}.}
    \item{Incoming Two-path (ITP)}{vertex \eqn{k} is an ITP shared partner of ordered pair \eqn{(i,j)} iff \eqn{j \to k \to i}{j->k->i}.}
    \item{Recursive Two-path (RTP)}{vertex \eqn{k} is a RTP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \leftrightarrow k \leftrightarrow j}{i<->k<->j}.}
    \item{Outgoing Shared Partner (OSP)}{vertex \eqn{k} is an OSP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k, j \to k}{i->k, j->k}.}
    \item{Incoming Shared Partner (ISP)}{vertex \eqn{k} is an ISP shared partner of ordered pair \eqn{(i,j)} iff \eqn{k \to i, k \to j}{k->i, k->j}.}
  }

Note that Robins et al. (2009) define closely related statistics to several of the above, using slightly different terminology.
}
%\value{\code{\link{ergm}} returns an object of class 'ergm' that is a list.
%}
\references{
Butts, CT.  (2008).  \dQuote{A Relational Event Framework for Social Action.} \emph{Sociological Methodology,} 38(1).

Hunter, D. R. and M. S. Handcock (2006). \dQuote{Inference in Curved Exponential Family Models for Networks.} \emph{Journal of Computational and Graphical Statistics,} 15: 565-583.

Robins, G; Pattison, P; and Wang, P.  (2009).  \dQuote{Closure, Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models for Directed Social Networks.}  \emph{Social Networks,} 31:105-117.

}
\seealso{statnet, network, ergm, ergm-terms}
\examples{
\donttest{
data(emon)
summary(emon[[1]]~desp(0:10,type="OSP"))
set.seed(1)
fit <- ergm(emon[[1]]~edges+mutual+dgwesp(0.2,type="OSP"), control=control.ergm(MCMC.interval=5000))
summary(fit)
}
}
\keyword{models}
