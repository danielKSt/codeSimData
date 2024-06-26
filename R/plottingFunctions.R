#' Plotting functions

#' Function for plotting a snapshot of some function/field on the surface
#' @param sField Matrix with surface field data for the interesting snapshot
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_tile theme
#' @export
surface.field.plot <- function(sField){
  data.loc <- expand.grid(X = 1:dim(sField)[1], Y = 1:dim(sField)[2])
  data.loc$sField <- array(sField)
  res <- ggplot(data=data.loc, mapping = aes(x = .data$X, y = .data$Y, fill = .data$sField))+
    geom_tile()+
    theme(legend.position = "bottom")
  return(res)
}


#' Function for plotting a snapshot of some point process on the surface
#' @param sPPCoords Coordinates for the point process
#' @importFrom ggplot2 ggplot aes geom_point
#' @export
surfacePP.plot <- function(sPPCoords){
  res <- ggplot(data = sPPCoords, mapping = aes(x = .data$x, y = .data$y))+
    geom_point()
  return(res)
}


#' Function for plotting a snapshot of surface field and surface pp together
#' @param sField Matrix with surface field data for the interesting snapshot
#' @param sPPCoords Coordinates for the point process
#'
#' @importFrom ggplot2 ggplot aes geom_tile theme geom_point
#'
#' @export
surface.field.pp.plot <- function(sField, sPPCoords){
  data.loc <- expand.grid(X = (1:dim(sField)[1])/256.0, Y = (1:dim(sField)[2])/256.0)
  data.loc$sField <- array(sField)
  res <- ggplot(data=data.loc, mapping = aes(x = .data$X, y = .data$Y, fill = .data$sField))+
    geom_tile()+
    theme(legend.position = "bottom")+
    geom_point(data = sPPCoords, mapping = aes(x = .data$x, y = .data$y, fill = .data$noPixels))
  return(res)
}
