InitErgmTerm.gwdspattr <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist,
                       varnames = c("attrname","p1keep","p2keep","alpha"),
                       vartypes = c("character","numeric", "numeric","numeric"),
                       defaultvalues = list(NULL, NULL, NULL, NULL),
                       required = c(TRUE, TRUE, TRUE, TRUE))  
  n <- network.size(nw)
  
  nodecov <- get.node.attr(nw, a$attrname)
  p1u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){ p1u<-c(p1u,NA) }
  nodecov <- match(nodecov,p1u,nomatch=length(p1u)+1)
  
  u <- c(a$p1keep,a$p2keep, a$alpha)
  
  coef.names <- paste0("gwdspattr.", a$attrname, ".", p1u[a$p1keep], "_", p1u[a$p2keep],".a:",a$alpha)
  list(name = "gwdspattr", coef.names = coef.names, #name and coef.names: required
       inputs = c(nodecov, u), minval=0)
}
