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
#' @param interp Set to TRUE if linear interpolation needs to be carried out
#' @param dr Interpolation parameter for radius
#' @param dh interpolation parameter for lag
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_contour_filled
#'
#' @export
aicPlot <- function(aicMatrix, lagValues, timeScale, lengthScale, radiar, interp = FALSE, dr = 1, dh = 1){
  if(interp){
    lMat_out <- akima::bilinear.grid(x = radiar, y = lagValues, z = aicMatrix, dx = dr, dy = dh)
    data.loc <- expand.grid(r = lMat_out$x*lengthScale, h = lMat_out$y*timeScale)
    data.loc$aic <- array(data = lMat_out$z)
  } else {
    data.loc <- expand.grid(r = radiar*lengthScale, h = lagValues*timeScale)
    data.loc$aic <- array(aicMatrix)
  }
  res <- ggplot(data=data.loc, mapping = aes(x = .data$r, y = .data$h, fill = .data$aic)) +
    geom_tile()
  return(res)
}


# setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData")
# library(codeSimData)
# library(ggplot2)
# library(viridis)
# library(akima)
#
# re <- 1000
# we <- "inf"
#
# settings <- paste("RE", re, "_WE", we, sep = "")
# load(file = paste("/Volumes/work/danielks/", settings, "/simulatedPoints.RDa", sep = ""))
# load(file = paste("/Volumes/work/danielks/", settings, "/scales.RDa", sep = ""))
# load(file = paste("/Volumes/work/danielks/", settings, "/res.RDa", sep = ""))
# res <- res_div[[1]]
#
# dimples_lags <- c(-4.0, 0.5)
# scars_lags <- c(-0.5, 0.5)
#
# lags <- res$lags
# radiar <- res$radiar
# dimples_ind.lags <- which((lags/timescale >= dimples_lags[1])*(lags/timescale <= dimples_lags[2]) == 1)
# scars_ind.lags <- which((lags/timescale >= scars_lags[1])*(lags/timescale <= scars_lags[2]) == 1)
#
# dimples <- res$dimples[dimples_ind.lags, ]
# dimpleOptim <- which(dimples[ , ] == max(dimples[ , ]), arr.ind = TRUE)
# temp_lags <- lags[dimples_ind.lags]
#
# lMat_out <- akima::bilinear.grid(x = radiar, y = temp_lags, z = t(dimples), dx = 1, dy = 1)
# data.loc <- expand.grid(r = lMat_out$x/timescale, h = lMat_out$y/lengthscale)
# data.loc$aic <- array(data = lMat_out$z)
#
# ggplot(data=data.loc, mapping = aes(x = .data$r, y = .data$h, fill = .data$aic)) +
#   geom_tile()
#
#
# dimplePlot <- aicPlot(aicMatrix = (t(dimples) - max(dimples)), lagValues = lags[dimples_ind.lags], timeScale = 1/timescale,
#                       lengthScale = 1/lengthscale, radiar = radiar, interp = TRUE) +
#   scale_fill_continuous(type = "viridis", direction = -1, name = bquote(Delta * "log(" * hat(L) * ")")) +
#   theme(legend.position = "right") +
#   geom_contour(mapping = aes(z = aic), color = "white", breaks = -c(min(res$global.dimples) - min(dimples)), linewidth = 1.2) +
#   geom_contour(mapping = aes(z = aic), bins = 5) +
#   geom_point(aes(x = radiar[dimpleOptim[2]]/lengthscale, y = temp_lags[dimpleOptim[1]]/timescale), color = "white") +
#   xlab("r") + ylab("h")
#
# dimplePlot
