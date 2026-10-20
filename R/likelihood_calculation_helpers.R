


#' Function for calculating the Poisson likelihood directly as a function of b
#' @description
#' Calculates the simultaneous likelihood for
#'
#' @param b Parameter value
#' @param logZ TxMxN array with the logarithm of the covariate
#' @param sumLogZ The sum of the array at the dimple/scar locations
#' @param n Complete count of points over all timesteps
#' @param pixelArea Area of a single pixel
#'
#' @export
calc_likelihood_poisson_b <- function(b, logZ, sumLogZ, n, pixelArea){
  # log(sum(Z^b)) is calculated from log(Z), since Z^b overflows or underflows when b is large
  return(b*sumLogZ-n*(log_sum_exp(b*logZ) + log(pixelArea)))
}

#' Function for calculating the complete Poisson likelihood
#' @description
#' Calculates the simultaneous likelihood for
#'
#' @param logA Logarithm of the parameter a, also named beta_0
#' @param b Parameter value - also named beta_1
#' @param logZ TxMxN array with the logarithm of the covariate
#' @param sumLogZ The sum of the array at the dimple/scar locations
#' @param n Complete count of points over all timesteps
#' @param pixelArea Area of a single pixel
#'
#' @export
calc_likelihood_poisson_complete <- function(logA, b, logZ, sumLogZ, n, pixelArea){
  return(n*logA + b*sumLogZ-exp(logA + log_sum_exp(b*logZ) + log(pixelArea)))
}

#' Find the maximum likelihood estimate of b
#' @description
#' Maximises the profile likelihood calc_likelihood_poisson_b. The profile likelihood is concave in b, so a maximum on
#' the edge of the search interval means that the estimate lies outside it, and the interval is then widened.
#'
#' @param logZ Array or vector with the logarithm of the covariate, one value per pixel
#' @param sumLogZ The sum of logZ at the dimple/scar locations
#' @param n Complete count of points over all timesteps
#' @param pixelArea Area covered by each value in logZ
#' @param bInterval Initial interval in which to search for the estimate
find_b_mle <- function(logZ, sumLogZ, n, pixelArea, bInterval){
  repeat{
    b <- stats::optimize(f = calc_likelihood_poisson_b, interval = bInterval, maximum = TRUE, pixelArea = pixelArea,
                         logZ = logZ, sumLogZ = sumLogZ, n = n)$maximum
    if(min(abs(b - bInterval)) > 1e-3*diff(bInterval)){
      break
    }
    if(max(abs(bInterval)) >= 1000){
      warning("The estimate of b is on the edge of the search interval (", bInterval[1], ", ", bInterval[2], ")")
      break
    }
    bInterval <- mean(bInterval) + 2*(bInterval - mean(bInterval))
  }
  return(b)
}

#' Calculate log(sum(exp(x))) without overflow or underflow
#'
#' @param x Numeric vector or array
log_sum_exp <- function(x){
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}


#' Function for calculating the sum of the covariate at the scar/dimple centroids
#' @description
#' Sums the values Z(x,y,t) of the covariate Z
#'
#' @param Z TxMxN array
#' @param structures List of data frames with the x and y coordinates of the scars/dimples, one data frame per
#' timestep in Z. Timesteps without any structures should have a data frame with zero rows.
#' @param l_x Length in x direction of window in which the structures are contained
#' @param l_y Length in y direction of window in which the structures are contained
#'
#' @export
get_sumLogZ_and_n <- function(Z, structures, l_x, l_y){
  sumLogZ <- 0
  n <- 0
  for (t in seq_along(structures)) {
    n <- n + nrow(structures[[t]])
    for(i in seq_len(nrow(structures[[t]]))){
      Zval <- get_Z_at_coord(Z, c(structures[[t]]$x[i], structures[[t]]$y[i]), t, l_x, l_y)
      sumLogZ <- sumLogZ + log(Zval)
    }
  }
  return(list(n = n, sumLogZ = sumLogZ))
}


#' Function for evaluating Z at a coordinate
#' @description
#' The array of covariates is not a continuous function, so we use the value of the pixel that contains the
#' coordinate. This matches the likelihood, where the integral over the window treats Z as constant within each pixel.
#' Pixel k in a dimension with N pixels covers ((k-1)*l/N, k*l/N).
#'
#' @param Z TxMxN array with covariate values
#' @param coord vector with (x,y)-coordinates
#' @param t Time at which Z should be evaluated
#' @param l_x Length in x direction of window in which the structures are contained
#' @param l_y Length in y direction of window in which the structures are contained
get_Z_at_coord <- function(Z, coord, t, l_x, l_y){
  if(coord[1] < 0 || coord[1] > l_x || coord[2] < 0 || coord[2] > l_y){
    stop("Coordinate (", coord[1], ", ", coord[2], ") is outside the observation window")
  }
  # A coordinate on the far edge of the window belongs to the last pixel
  xInd <- min(floor(coord[1]*dim(Z)[2]/l_x) + 1, dim(Z)[2])
  yInd <- min(floor(coord[2]*dim(Z)[3]/l_y) + 1, dim(Z)[3])
  return(Z[t, xInd, yInd])
}
