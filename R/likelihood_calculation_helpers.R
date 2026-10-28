


#' Function for calculating the Poisson likelihood directly as a function of b
#' @description
#' Calculates the simultaneous likelihood for the model intensity = a*Z^b, maximised over a. With Z in place of logZ it
#' is the likelihood of the model with the intensity on log scale, log(intensity) = beta_0 + b*Z
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
#' Calculates the simultaneous likelihood for the model intensity = a*Z^b. With Z in place of logZ it is the likelihood
#' of the model with the intensity on log scale, log(intensity) = logA + b*Z
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


#' Function for evaluating Z at all the scar/dimple centroids
#' @description
#' The values Z(x,y,t) at each scar/dimple, with the same pixel convention as get_sumLogZ_and_n
#'
#' @param Z TxMxN array
#' @param structures List of data frames with the x and y coordinates of the scars/dimples, one data frame per
#' timestep in Z
#' @param l_x Length in x direction of window in which the structures are contained
#' @param l_y Length in y direction of window in which the structures are contained
#'
#' @returns Vector with one value per point, ordered by timestep
get_Z_at_points <- function(Z, structures, l_x, l_y){
  # as.numeric gives numeric(0) instead of NULL when there are no points
  return(as.numeric(unlist(lapply(seq_along(structures), function(t){
    vapply(seq_len(nrow(structures[[t]])), function(i){
      get_Z_at_coord(Z, c(structures[[t]]$x[i], structures[[t]]$y[i]), t, l_x, l_y)
    }, numeric(1))
  }))))
}


# ----------------------------------------------------------------------------------------------------------------------
# Added baseline intensity
#
# The model nu(x,y,t) = nu_0 + a*Z(x,y,t)^b = a*(lambda + Z(x,y,t)^b), with b >= 0. Z may be zero, and Z^b is then zero
# also for b = 0, which is the limit from b > 0. With A the total area of all snapshots and S_b the integral of Z^b over
# all snapshots, the maximum likelihood estimates satisfy nu_0*A + a*S_b = n, the number of points. Instead of lambda,
# the functions below use the share w = nu_0*A/n of the expected count that comes from the baseline, so that
#   nu_0 = w*n/A,  a = (1-w)*n/S_b,  lambda = w*S_b/((1-w)*A).
# The profile log likelihood is then
#   l(w, b) = n*log(n) - n + sum_i log(w/A + (1-w)*Z_i^b/S_b),
# which is concave in w for fixed b. w lies in [0, 1], where w = 0 is the model without baseline and w = 1 a constant
# intensity.
# ----------------------------------------------------------------------------------------------------------------------

#' Function for calculating the complete Poisson likelihood with an added baseline intensity
#' @description
#' Calculates the log likelihood of the model nu = nu_0 + a*Z^b, up to the same constant as
#' calc_likelihood_poisson_complete
#'
#' @param nu0 Baseline intensity
#' @param a Parameter a
#' @param b Parameter b, at least 0
#' @param logZ Logarithm of the covariate at the pixels where it is positive
#' @param logZAtPoints Logarithm of the covariate at the dimple/scar locations, -Inf where the covariate is zero
#' @param pixelArea Area of a single pixel
#' @param totalArea Area of all the pixels, also those where the covariate is zero
#'
#' @export
calc_likelihood_baseline_complete <- function(nu0, a, b, logZ, logZAtPoints, pixelArea, totalArea){
  # log(nu_0 + a*Z^b) and the integral a*S_b are calculated from log(Z), since Z^b overflows or underflows when b is
  # large
  logNu <- log_add_exp(log(nu0), log(a) + power_log(b, logZAtPoints))
  return(sum(logNu) - nu0*totalArea - exp(log(a) + log(pixelArea) + log_sum_exp(b*logZ)))
}

#' Function for calculating the profile Poisson likelihood with an added baseline intensity as a function of b
#' @description
#' Calculates the log likelihood for a given b, maximised over nu_0 and a, see find_w_mle
#'
#' @param b Parameter value, at least 0
#' @param logZ Logarithm of the covariate at the pixels where it is positive
#' @param logZAtPoints Logarithm of the covariate at the dimple/scar locations, -Inf where the covariate is zero
#' @param pixelArea Area of a single pixel
#' @param totalArea Area of all the pixels, also those where the covariate is zero
#'
#' @export
calc_likelihood_baseline_b <- function(b, logZ, logZAtPoints, pixelArea, totalArea){
  return(find_w_mle(b = b, logZ = logZ, logZAtPoints = logZAtPoints, pixelArea = pixelArea,
                    totalArea = totalArea)$logLik)
}

#' Find the maximum likelihood estimate of w for a given b
#' @description
#' The profile likelihood is concave in w, so the estimate is 0 or 1 when the derivative at that end of the interval
#' from 0 to 1 points out of the interval, and otherwise the maximum inside the interval
#'
#' @param b Parameter value, at least 0
#' @param logZ Logarithm of the covariate at the pixels where it is positive
#' @param logZAtPoints Logarithm of the covariate at the dimple/scar locations, -Inf where the covariate is zero
#' @param pixelArea Area of a single pixel
#' @param totalArea Area of all the pixels, also those where the covariate is zero
#'
#' @returns List with w, log(S_b) and the log likelihood
find_w_mle <- function(b, logZ, logZAtPoints, pixelArea, totalArea){
  n <- length(logZAtPoints)
  logS <- log(pixelArea) + log_sum_exp(b*logZ)
  if(n == 0){
    # Without points the estimates of nu_0 and a are 0, and so is the log likelihood
    return(list(w = NA_real_, logS = logS, logLik = 0))
  }
  # With c_i = Z_i^b/S_b and alpha = 1/A the profile likelihood is n*log(n) - n + sum(log(w*alpha + (1-w)*c_i))
  logC <- power_log(b, logZAtPoints) - logS
  alpha <- 1/totalArea
  profile <- function(w){
    if(w == 0) sum(logC) else sum(log(w*alpha + (1-w)*exp(logC)))
  }
  # The derivative is alpha*sum(1/c_i) - n at w = 0, which is infinite when a point is on a zero pixel, and
  # n - sum(c_i)/alpha at w = 1
  if(alpha*sum(exp(-logC)) <= n){
    w <- 0
  } else if(sum(exp(logC))/alpha <= n){
    w <- 1
  } else {
    w <- stats::optimize(f = profile, interval = c(0, 1), maximum = TRUE, tol = 1e-10)$maximum
  }
  return(list(w = w, logS = logS, logLik = n*log(n) - n + profile(w)))
}

#' Find the maximum likelihood estimates of b and w
#' @description
#' The profile likelihood in b does not have to be concave, so it is first evaluated at bGrid, and then maximised
#' between the neighbours of the best grid value. The grid is extended when the best value is its largest
#'
#' @param logZ Logarithm of the covariate at the pixels where it is positive
#' @param logZAtPoints Logarithm of the covariate at the dimple/scar locations, -Inf where the covariate is zero
#' @param pixelArea Area of a single pixel
#' @param totalArea Area of all the pixels, also those where the covariate is zero
#' @param bGrid Values of b, at least 0, at which the profile likelihood is evaluated first
#'
#' @returns List with b, w, log(S_b) and the log likelihood
find_baseline_mle <- function(logZ, logZAtPoints, pixelArea, totalArea, bGrid){
  if(any(bGrid < 0)){
    stop("b must be at least 0 in the model with a baseline intensity")
  }
  profile <- function(b){
    calc_likelihood_baseline_b(b = b, logZ = logZ, logZAtPoints = logZAtPoints, pixelArea = pixelArea,
                               totalArea = totalArea)
  }
  bGrid <- sort(unique(bGrid))
  values <- vapply(bGrid, profile, numeric(1))
  while(which.max(values) == length(bGrid)){
    if(max(bGrid) >= 1000){
      warning("The estimate of b is at the largest value tried, ", max(bGrid))
      break
    }
    bGrid <- c(bGrid, 2*max(bGrid, 0.5))
    values <- c(values, profile(max(bGrid)))
  }
  k <- which.max(values)
  opt <- stats::optimize(f = profile, interval = bGrid[c(max(k - 1, 1), min(k + 1, length(bGrid)))], maximum = TRUE)
  # optimize does not evaluate the ends of the interval, so the best grid value is kept if it is better
  b <- if(opt$objective > values[k]) opt$maximum else bGrid[k]
  return(c(list(b = b), find_w_mle(b = b, logZ = logZ, logZAtPoints = logZAtPoints, pixelArea = pixelArea,
                                   totalArea = totalArea)))
}

#' Calculate b*log(Z) where Z^b is zero for Z = 0, also for b = 0
#'
#' @param b Parameter value, at least 0
#' @param logZ Numeric vector with the logarithm of the covariate, -Inf where it is zero
power_log <- function(b, logZ){
  return(ifelse(logZ == -Inf, -Inf, b*logZ))
}

#' Calculate log(exp(x) + exp(y)) without overflow or underflow
#'
#' @param x Numeric vector
#' @param y Numeric vector
log_add_exp <- function(x, y){
  m <- pmax(x, y)
  res <- m + log1p(exp(-abs(x - y)))
  # Where both are -Inf, x - y is NaN
  res[m == -Inf] <- -Inf
  return(res)
}


# ----------------------------------------------------------------------------------------------------------------------
# Intensity on log scale
#
# The standard model of Poisson regression, log(nu(x,y,t)) = beta_0 + beta_1*Z(x,y,t). Unlike the power model
# nu = a*Z^b = exp(beta_0 + b*log(Z)), Z may be zero or negative. The power model is this model with log(Z) in place of
# Z, so calc_likelihood_poisson_b and calc_likelihood_poisson_complete give the log likelihood with Z in place of logZ,
# with the profile likelihood
#   l(beta_1) = beta_1*sum_i Z_i - n*log(S(beta_1)) + n*log(n) - n,  S(beta_1) = integral of exp(beta_1*Z),
# and exp(beta_0) = n/S(beta_1). The profile likelihood is concave in beta_1, so it has a single maximum.
# ----------------------------------------------------------------------------------------------------------------------

#' Find the maximum likelihood estimate of beta_1 in the model with the intensity on log scale
#' @description
#' Maximises the profile likelihood calc_likelihood_poisson_b with Z in place of logZ. The scale of beta_1 is 1/Z, so
#' the search is done for beta_1*sd(Z), which does not depend on the units of Z
#'
#' @param Z Array or vector with the covariate, one value per pixel
#' @param sumZ The sum of Z at the dimple/scar locations
#' @param n Complete count of points over all timesteps
#' @param pixelArea Area covered by each value in Z
#' @param betaInterval Initial interval in which to search for beta_1*sd(Z), see find_b_mle
#'
#' @returns The estimate of beta_1. It is 0 when Z is constant, where beta_1 can not be estimated
find_beta1_mle <- function(Z, sumZ, n, pixelArea, betaInterval){
  scaleZ <- stats::sd(as.vector(Z))
  if(!(scaleZ > 0)){
    return(0)
  }
  beta1Scaled <- find_b_mle(logZ = Z/scaleZ, sumLogZ = sumZ/scaleZ, n = n, pixelArea = pixelArea,
                            bInterval = betaInterval)
  return(beta1Scaled/scaleZ)
}
