.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.terms.contrib", c("statnet"), FALSE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}



