
#' @useDynLib codeSimData distanceOnTorus_

distanceOnTorus <- function(d, x, y){
  distance <- .C("distanceOnTorus_", as.integer(d), as.double(x), as.double(y), result = as.double(0))
  return(as.double(distance$result))
}
