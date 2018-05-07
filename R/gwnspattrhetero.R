InitErgmTerm.gwnspattrhetero <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=TRUE, 
                       varnames = c("attrname","pexcl","alpha"),
                       vartypes = c("character","numeric","numeric"),
                       defaultvalues = list(NULL, NULL, NULL),
                       required = c(TRUE, TRUE, TRUE))  
  n <- network.size(nw)
  
  nodecov <- get.node.attr(nw, a$attrname)
  exclu<-sort(unique(nodecov))
  if(any(is.na(nodecov))){ exclu<-c(exclu,NA) }
  nodecov <- match(nodecov,exclu,nomatch=length(exclu)+1)
  
  u <- c(a$pexcl, a$alpha)
  
  coef.names <- paste0("gwnspattrhetero.", a$attrname, ".excl:", exclu[a$pexcl], ".a:",a$alpha)
  list(name = "gwnspattrhetero", coef.names = coef.names, #name and coef.names: required
       inputs = c(nodecov, u), minval=0)
}
