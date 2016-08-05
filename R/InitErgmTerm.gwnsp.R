#Term to count NSP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Recursive two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per NSP term.  The default, OTP, retains
#the original behavior of nsp/gwnsp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#
InitErgmTerm.dnsp<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type"),
                      vartypes = c("numeric","character"),
                      defaultvalues = list(NULL,"OTP"),
                      required = c(TRUE, FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for sp; valid types are:",paste(type.vec, collapse=","))
  dname<-"nsp"
  if(is.directed(nw)){
    conam <- paste("nsp",type,sep=".")
    typecode<-which(type==type.vec)
    dname <- "dnsp"
  }else{
    message("Use the ergm term 'nsp' for undirected networks.")
    conam<-"nsp"
    type<-"UTP"
    typecode<-0
  }
  list(name=dname, coef.names=paste(conam,d,sep=""), inputs=c(typecode,d), minval=0)
}


################################################################################
InitErgmTerm.dgwnsp<-function(nw, arglist, initialfit=FALSE, ...) {
  # the following line was commented out in <InitErgm.gwnsp>:
  #    ergm.checkdirected("gwnsp", is.directed(nw), requirement=FALSE)
  # so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("alpha","fixed","cutoff","type"),
                      vartypes = c("numeric","logical","numeric","character"),
                      defaultvalues = list(0, FALSE, 30,"OTP"),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  alpha<-a$alpha;fixed<-a$fixed
  cutoff<-a$cutoff
  alpha=alpha[1] # Not sure why anyone would enter a vector here, but...
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code; valid types are:",paste(type.vec, collapse=","))
  dname<-"dnsp"
  
  if(!is.directed(nw)){  
    stop("Use the gwnsp term for undirected networks.")
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwnsp",type,sep=".")
  }
  
  if(!initialfit && !fixed){ # This is a curved exponential family model
    #   d <- 1:(network.size(nw)-1)
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...){
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...){
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    
    params<-list(gwnsp=NULL,gwnsp.alpha=alpha)
    names(params)<-c(basenam,paste(basenam,"alpha",sep="."))
    
    list(name=dname, coef.names=paste("nsp#",d,sep=""),
         inputs=c(typecode,d), params=params,
         map=map, gradient=gradient)
  }else{
    dname<-"dgwnsp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if (initialfit && !fixed)  # First pass to get MPLE coefficient
      coef.names <- basenam
    else { # fixed == TRUE
      if (is.directed(nw)) 
        coef.names <- paste("gwnsp",type,"fixed",alpha,sep=".")
      else
        coef.names <- paste("gwnsp.fixed",alpha,sep=".")
    }
    
    list(name=dname, coef.names=coef.names, inputs=c(alpha,typecode,maxesp))
  }
}
  