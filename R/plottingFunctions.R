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
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_contour_filled
#'
#' @export
aicPlot <- function(aicMatrix, lagValues, timeScale, lengthScale, radiar){
  filled <- fill_missing(mat_in = aicMatrix, x_in = radiar, y_in = lagValues)
  aicMatrix <- filled$mat
  radiar <- filled$x
  lagValues <- filled$y
  data.loc <- expand.grid(r = radiar*lengthScale, h = lagValues*timeScale)
  data.loc$aic <- array(aicMatrix)
  res <- ggplot(data=data.loc, mapping = aes(x = .data$r, y = .data$h, fill = .data$aic)) +
    geom_tile()
  return(res)
}


#' Function for filling in missing values in a image matrix using linear regression
#' This assumes all the different dx/dy values are multiples of eachother
#' @param mat_in Matrix as input
#' @param x_in Values connected to the x-indices
#' @param y_in Values connected to the y-indices
#'
fill_missing <- function(mat_in, x_in, y_in){
  dx_max <- x_in[2] - x_in[1]
  dx_min <- x_in[2] - x_in[1]
  for (i in 3:length(x_in)) {
    dx <- x_in[i] - x_in[i-1]
    dx_max <- max(dx, dx_max)
    dx_min <- min(dx, dx_min)
  }
  if(dx_max > dx_min){
    x_out <- seq(from = min(x_in), to = max(x_in), by = dx_min)
    mat_out <- matrix(data = 0, ncol = length(y_in), nrow = length(x_out))
    mat_out[1, ] <- mat_in[1, ]
    i_out <- 2
    for (i in 2:length(x_in)) {
      if(x_in[i]- x_in[i-1] > dx_min){
        n <- as.integer((x_in[i] - x_in[i-1])/dx_min)
        for(k in 1:n){
          mat_out[i_out, ] <- (k*mat_in[i, ] + (n-k)*mat_in[i-1, ])/n
          i_out <- i_out + 1
        }
      } else {
        mat_out[i_out, ] <- mat_in[i, ]
        i_out <- i_out + 1
      }
    }
    mat_in <- mat_out
  } else {
    x_out <- x_in
  }

  dy_max <- y_in[2] - y_in[1]
  dy_min <- y_in[2] - y_in[1]
  for (i in 3:length(y_in)) {
    dy <- y_in[i] - y_in[i-1]
    dy_max <- max(dy, dy_max)
    dy_min <- min(dy, dy_min)
  }
  if(dy_max > dy_min){
    y_out <- seq(from = min(y_in), to = max(y_in), by = dy_min)
    mat_out <- matrix(data = 0, ncol = length(y_in), nrow = length(y_out))
    mat_out[ , 1] <- mat_in[ , 1]
    i_out <- 2
    for (i in 2:length(y_in)) {
      if(y_in[i]- y_in[i-1] > dy_min){
        n <- as.integer((y_in[i] - y_in[i-1])/dy_min)
        for(k in 1:n){
          mat_out[ , i_out] <- (k*mat_in[ ,i] + (n-k)*mat_in[ , i-1])/n
          i_out <- i_out + 1
        }
      } else {
        mat_out[ , i_out] <- mat_in[ , i]
        i_out <- i_out + 1
      }
    }
  } else {
    y_out <- y_in
    mat_out <- mat_in
  }

  return(list(mat = mat_out, x = x_out, y = y_out))
}

#
# library(ggplot2)
# x_in <- c(1:5, seq(from = 6, to = 17, by = 3))
# y_in <- c(1:10)
#
# mat_in <- matrix(data = 0, nrow = length(x_in), ncol = length(y_in))
# for (x_ind in 1:length(x_in)) {
#   for(y_ind in y_in){
#     mat_in[x_ind, y_ind] <- x_in[x_ind]*y_ind
#   }
# }
#
#
# data.loc <- expand.grid(r = x_in, h = y_in)
# data.loc$fill <- array(mat_in)
# ggplot(data=data.loc, mapping = aes(x = r, y = h, fill = fill)) + geom_tile()
#
# a <- aicPlot(aicMatrix = mat_in, lagValues = y_in, radiar = x_in, lengthScale = 1, timeScale = 1)
# a

