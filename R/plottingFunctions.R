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
  data.loc <- expand.grid(X = (1:dim(sField)[1]), Y = (1:dim(sField)[2]))
  data.loc$sField <- array(sField)
  res <- ggplot(data=data.loc, mapping = aes(x = .data$X, y = .data$Y, fill = .data$sField))+
    geom_tile()+
    theme(legend.position = "bottom")+
    geom_point(data = sPPCoords, mapping = aes(x = .data$x*256, y = .data$y*256, fill = .data$noPixels))
  return(res)
}


#' Function for plotting a snapshot of surface field and surface pp together
#' @param aicMatrix Matrix with surface field data for the interesting snapshot
#' @param lagValues Coordinates for the point process
#' @param timeScale how many timescales per frame?
#' @param lengthScale how many lengthscales per pixle?
#' @param radiar description
#' @param plotAsContour set to TRUE if you want a contour plot instead of a heatmap
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_contour_filled
#'
#' @export
aicPlot <- function(aicMatrix, lagValues, timeScale, lengthScale, radiar, plotAsContour = FALSE){
  data.loc <- expand.grid(r = radiar*lengthScale, h = lagValues*timeScale)
  data.loc$aic <- array(aicMatrix)
  res <- ggplot(data=data.loc, mapping = aes(x = .data$r, y = .data$h, fill = .data$aic, z = .data$aic))
  if(plotAsContour){
    res <- res + geom_contour_filled()
  }
  else {
    res <- res + geom_tile()
  }
  return(res)
}
