#' Functions used for manipulation of the surface divergence


#' Function for estimating the norm of the gradient field of the surface divergence
#'
#' @param sDiv surface divergence
#' @param xdim Dimension in the x direction
#' @param ydim Dimension in the y direction
#' @param hx Step size in x-direction
#' @param hy Step size in y-direction
#'
#' @returns Matrix with the norm of the gradient of the surface divergence
#'
#' @export
estimate.gradient.norm.sdiv <- function(sDiv, xdim, ydim, hx, hy){
  res <- matrix(data = NA, nrow = xdim, ncol = ydim)
  # Corner pixels
  res[1,1] <- sqrt(((sDiv[2,1]-sDiv[xdim,1])/(2*hx))^2+((sDiv[1, 2]-sDiv[1,ydim])/(2*hy))^2)
  res[xdim,1] <- sqrt(((sDiv[1,1]-sDiv[xdim-1,1])/(2*hx))^2+((sDiv[xdim, 2]-sDiv[xdim,ydim])/(2*hy))^2)
  res[1,ydim] <- sqrt(((sDiv[2,ydim]-sDiv[xdim,ydim])/(2*hx))^2+((sDiv[1, 1]-sDiv[1,ydim-1])/(2*hy))^2)
  res[xdim,ydim] <- sqrt(((sDiv[1,ydim]-sDiv[xdim-1,ydim])/(2*hx))^2+((sDiv[xdim, 1]-sDiv[xdim,ydim-1])/(2*hy))^2)

  # Edge rows
  for (k in 2:(xdim-1)) {
    res[k, 1] <- sqrt(((sDiv[k+1,1]-sDiv[k-1,1])/(2*hx))^2+((sDiv[k,2]-sDiv[k,ydim])/(2*hy))^2)
    res[k, ydim] <- sqrt(((sDiv[k+1,ydim]-sDiv[k-1,ydim])/(2*hx))^2+((sDiv[k,1]-sDiv[k,ydim-1])/(2*hy))^2)
  }

  for (k in 2:(ydim-1)) {
    res[1, k] <- sqrt(((sDiv[2, k]-sDiv[xdim,k])/(2*hx))^2+((sDiv[1,k+1]-sDiv[1,k-1])/(2*hy))^2)
    res[xdim, k] <- sqrt(((sDiv[1, k]-sDiv[xdim-1,k])/(2*hx))^2+((sDiv[xdim,k+1]-sDiv[xdim,k-1])/(2*hy))^2)
  }

  # Central pixels
  for (i in 2:(xdim-1)) {
    for (j in 2:(ydim-1)) {
      res[i,j] <- sqrt(((sDiv[i+1,j]-sDiv[i-1,j])/(2*hx))^2+
                         ((sDiv[i,j+1]-sDiv[i,j-1])/(2*hy))^2)
    }
  }
  return(res)
}


#' Funciton for smoothing out a single datapoint
#' @param sField Matrix of values on the surface
#' @param vindauga size of slinding window for smoothing
#' @param glatter matrix of weights for how to smooth the pixels
#' @param i First index
#' @param j Second index
#'
#' @returns A smoothed value
single.pixle.smoothing<- function(sField, vindauga, glatter, i, j){
  smoother.input <- c(1:dim(vindauga)[1])
  for (k in c(1:dim(vindauga)[1])) {
    smoother.input[k] <- sField[(i-1-vindauga[k,1])%%dim(sField)[1]+1, (j-1-vindauga[k,2])%%dim(sField)[2]+1]
  }
  return(glatter(smoother.input))
}

#' Function for smoothing out data for some field on the surface
#'
#' @param sField Matrix of values on the surface
#' @param vindauga size of slinding window for smoothing
#' @param glatter matrix of weights for how to smooth the pixels
#'
#' @returns A smoothed version of sField
#'
#' @export
smooth.surface.field <- function(sField, vindauga, glatter){
  xdim <- dim(sField)[1]
  ydim <- dim(sField)[2]
  res <- matrix(data = NA, nrow = xdim, ncol = ydim)
  for (i in 1:xdim) {
    for (j in 1:ydim) {
      res[i,j] <- single.pixle.smoothing(sField, vindauga, glatter, i, j)
    }
  }
  return(res)
}
