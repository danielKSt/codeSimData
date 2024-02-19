#' Distance between two points on a flat torus of arbitrary dimension
#'
#' @param d The dimension of the torus
#' @param x The first point
#' @param y The second point
#'
#' @return The distance between x and y.
#' @export
#'
#' @examples
#' d <- 1
#' x <- 0.1
#' y <- 0.2
#' distanceOnTorus(d, x, y)
#'
#' @useDynLib codeSimData distanceOnTorus_

distanceOnTorus <- function(d, x, y){
  distance <- .C("distanceOnTorus_", as.integer(d), as.double(x), as.double(y), result = as.double(0))
  return(as.double(distance$result))
}
