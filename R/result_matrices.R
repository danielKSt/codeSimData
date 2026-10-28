
#' Prepare the points, given a list of all points
#' @description
#' Organize the points into a list with the desired time indices, scaled from the unit square to the correct dimensions
#'
#' @param pointsList List of all points in a unit square
#' @param l_x Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param l_y Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param tidspunkt Timesteps to use for list
#'
#' @export
make.points.list <- function(pointsList, l_x, l_y, tidspunkt){
  res <- vector(mode = "list", length = length(tidspunkt))
  if(is.numeric(tidspunkt)){
    for(i in 1:length(tidspunkt)){
      t <- tidspunkt[i]
      if(is.empty.snapshot(pointsList[[t]])){
        # Assigning NULL would delete the element, so empty timesteps get a data frame with zero rows
        res[[i]] <- data.frame(x = numeric(0), y = numeric(0))
      } else {
        res[[i]] <- data.frame(x = pointsList[[t]]$x*l_x,
                               y = pointsList[[t]]$y*l_y)
      }
    }
  } else {
    for(i in 1:length(tidspunkt)){
      ensemble <- tidspunkt[[i]][1]
      t <- tidspunkt[[i]][2]
      if(is.empty.snapshot(pointsList[[ensemble]][[t]])){
        res[[i]] <- data.frame(x = numeric(0), y = numeric(0))
      } else {
        res[[i]] <- data.frame(x = pointsList[[ensemble]][[t]]$x*l_x,
                               y = pointsList[[ensemble]][[t]]$y*l_y)
      }
    }
  }
  return(res)
}


#' Fit the regression model from the paper
#' @description
#' The paper describes a regression model given a certain covariate Z, the regression model is then
#' intensity(x,y,t)=a*Z(x,y,t)^b,
#' where Z might be the local variance, or some other local flow property.
#'
#' @param Z TxMxN array of covariates
#' @param structures List of scars/dimples locations, only at the time indices that we want to use
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param returnParams Set to TRUE to return the values a and b
#' @param bInterval Initial interval in which to search for the maximum likelihood estimate of b. The interval is
#' widened if the maximum is on its edge
#'
#' @export
fit.given.covar <- function(Z, structures, l_x, l_y, returnParams, bInterval = c(-20, 20)){
  # Every pixel of every snapshot is one term in sum(Z^b), so each term is weighted by the area of a single pixel
  pixelArea <- l_x*l_y/(dim(Z)[2]*dim(Z)[3])
  logZ <- log(Z)
  if(any(!is.finite(logZ))){
    stop("The covariate Z must be positive and finite, since the model uses log(Z)")
  }
  inData <- get_sumLogZ_and_n(Z = Z, structures = structures, l_x = l_x, l_y = l_y)
  b <- find_b_mle(logZ = logZ, sumLogZ = inData$sumLogZ, n = inData$n, pixelArea = pixelArea, bInterval = bInterval)
  # a = n/(sum(Z^b)*pixelArea), calculated on the log scale
  logA <- log(inData$n) - log(pixelArea) - log_sum_exp(b*logZ)
  logLik <- calc_likelihood_poisson_complete(logA, b, logZ, inData$sumLogZ, inData$n, pixelArea)
  if(returnParams){
    return(list(a = exp(logA), b = b, logLik = logLik))
  } else {
    return(logLik)
  }
}

#' Fit the regression model with an added baseline intensity
#' @description
#' The regression model with an added baseline intensity nu_0 is
#' intensity(x,y,t)=nu_0+a*Z(x,y,t)^b=a*(lambda+Z(x,y,t)^b),
#' with b >= 0. Unlike in fit.given.covar, Z may be zero, and Z^b is then zero, also for b = 0. With nu_0 = 0 this is
#' the model of fit.given.covar, so the log likelihood is at least as large as from fit.given.covar when its estimate of
#' b is at least 0.
#'
#' @param Z TxMxN array of non-negative covariates
#' @param structures List of scars/dimples locations, only at the time indices that we want to use
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param returnParams Set to TRUE to return the parameters
#' @param bGrid Values of b, at least 0, at which the profile likelihood is evaluated before it is maximised. The
#' profile likelihood does not have to be concave in b, so the grid should be fine enough to find the right maximum
#'
#' @returns The log likelihood, or a list with a, b, nu0, lambda = nu0/a, the share w of the expected number of points
#' that comes from the baseline, and the log likelihood. When b = 0 and Z has no zeros the intensity is constant for
#' every w, and w = 0 is returned
#'
#' @export
fit.given.covar.baseline <- function(Z, structures, l_x, l_y, returnParams,
                                     bGrid = c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8)){
  if(any(!is.finite(Z)) || any(Z < 0)){
    stop("The covariate Z must be non-negative and finite")
  }
  if(all(Z == 0)){
    stop("The covariate Z is zero everywhere")
  }
  pixelArea <- l_x*l_y/(dim(Z)[2]*dim(Z)[3])
  totalArea <- pixelArea*length(Z)
  logZAtPoints <- log(get_Z_at_points(Z = Z, structures = structures, l_x = l_x, l_y = l_y))
  fit <- find_baseline_mle(logZ = log(Z[Z > 0]), logZAtPoints = logZAtPoints, pixelArea = pixelArea,
                           totalArea = totalArea, bGrid = bGrid)
  if(!returnParams){
    return(fit$logLik)
  }
  n <- length(logZAtPoints)
  nu0 <- fit$w*n/totalArea
  a <- (1 - fit$w)*n/exp(fit$logS)
  return(list(a = a, b = fit$b, nu0 = nu0, lambda = nu0/a, w = fit$w, logLik = fit$logLik))
}

#' Fit the regression model with the intensity on log scale
#' @description
#' The standard model of Poisson regression,
#' log(intensity(x,y,t))=beta_0+beta_1*Z(x,y,t).
#' Unlike in fit.given.covar, Z may be zero or negative, for example the positive or negative part of the divergence,
#' or the divergence itself. The intensity is exponential in Z, while fit.given.covar gives an intensity that is a power
#' of Z, which is this model with log(Z) in place of Z. Both models have two parameters.
#'
#' @param Z TxMxN array of covariates
#' @param structures List of scars/dimples locations, only at the time indices that we want to use
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param returnParams Set to TRUE to return the values beta_0 and beta_1
#' @param betaInterval Initial interval in which to search for the maximum likelihood estimate of beta_1, in units of
#' 1/sd(Z). The interval is widened if the maximum is on its edge
#'
#' @returns The log likelihood, or a list with beta0, beta1 and the log likelihood. When Z is constant beta_1 can not be
#' estimated, and beta1 = 0 is returned
#'
#' @export
fit.given.covar.loglinear <- function(Z, structures, l_x, l_y, returnParams, betaInterval = c(-20, 20)){
  if(any(!is.finite(Z))){
    stop("The covariate Z must be finite")
  }
  pixelArea <- l_x*l_y/(dim(Z)[2]*dim(Z)[3])
  sumZ <- sum(get_Z_at_points(Z = Z, structures = structures, l_x = l_x, l_y = l_y))
  n <- sum(vapply(structures, nrow, integer(1)))
  beta1 <- find_beta1_mle(Z = Z, sumZ = sumZ, n = n, pixelArea = pixelArea, betaInterval = betaInterval)
  # exp(beta_0) = n/(sum(exp(beta_1*Z))*pixelArea), calculated on the log scale
  beta0 <- log(n) - log(pixelArea) - log_sum_exp(beta1*Z)
  logLik <- calc_likelihood_poisson_complete(logA = beta0, b = beta1, logZ = Z, sumLogZ = sumZ, n = n,
                                             pixelArea = pixelArea)
  if(returnParams){
    return(list(beta0 = beta0, beta1 = beta1, logLik = logLik))
  } else {
    return(logLik)
  }
}

#' Fit the global model from the paper
#' @description
#' In the global model the covariate is constant over each snapshot, Z(x,y,t)=Z_t, so the intensity is
#' intensity(x,y,t)=a*Z_t^b.
#' The likelihood then only depends on Z_t and the number of points in each snapshot. This gives the same result as
#' fit.given.covar with a covariate that is constant over each snapshot, but is much faster.
#'
#' @param Zglobal Vector with the covariate value of each snapshot
#' @param structures List of scars/dimples locations, one element per snapshot in Zglobal
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param returnParams Set to TRUE to return the values a and b
#' @param bInterval Initial interval in which to search for the maximum likelihood estimate of b. The interval is
#' widened if the maximum is on its edge
#'
#' @export
fit.global.covar <- function(Zglobal, structures, l_x, l_y, returnParams, bInterval = c(-20, 20)){
  if(length(Zglobal) != length(structures)){
    stop("Zglobal and structures must have one element per snapshot")
  }
  logZ <- log(Zglobal)
  if(any(!is.finite(logZ))){
    stop("The covariate Z must be positive and finite, since the model uses log(Z)")
  }
  counts <- vapply(structures, nrow, integer(1))
  n <- sum(counts)
  sumLogZ <- sum(counts*logZ)
  # Each snapshot acts as a single pixel that covers the whole window
  winArea <- l_x*l_y
  b <- find_b_mle(logZ = logZ, sumLogZ = sumLogZ, n = n, pixelArea = winArea, bInterval = bInterval)
  # a = n/(sum(Zglobal^b)*winArea), calculated on the log scale
  logA <- log(n) - log(winArea) - log_sum_exp(b*logZ)
  logLik <- calc_likelihood_poisson_complete(logA, b, logZ, sumLogZ, n, winArea)
  if(returnParams){
    return(list(a = exp(logA), b = b, logLik = logLik))
  } else {
    return(logLik)
  }
}

#' Collect the point sets that are given
#' @description
#' Puts the point sets that are not NULL in a named list, in the order scars, dimples, detectedDimples
#'
#' @param scars List of scars, or NULL
#' @param dimples List of dimples, or NULL
#' @param detectedDimples List of dimples detected from the surface elevation, or NULL
#'
#' @returns Named list with the point sets that are not NULL
collect.point.sets <- function(scars, dimples, detectedDimples){
  pointSets <- list(scars = scars, dimples = dimples, detectedDimples = detectedDimples)
  pointSets <- pointSets[!vapply(pointSets, is.null, logical(1))]
  if(length(pointSets) == 0){
    stop("At least one of scars, dimples and detectedDimples must be given")
  }
  return(pointSets)
}

#' Calculate the log likelihood of the global model for a specific lag
#' @description
#' Compute the log likelihood of the global model, where the covariate is a statistic of the whole snapshot
#'
#' @param h Lag, in timesteps, from the snapshots of the points to the snapshots of the field
#' @param field TxMxN array with the field the covariate is calculated from
#' @param scars List of scars at each timestep in tidspunkt, as returned by make.points.list, or NULL to leave them out
#' @param dimples List of dimples at each timestep in tidspunkt, as returned by make.points.list, or NULL to leave them
#' out
#' @param detectedDimples List of dimples detected from the surface elevation at each timestep in tidspunkt, as returned
#' by make.points.list, or NULL to leave them out
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param tidspunkt Vector with the timesteps of the points
#' @param localStat Statistic of the snapshot used as covariate. "meanAbs" gives the mean absolute value of the
#' snapshot, and "variance" and "meanSquare" both give the mean square
#'
#' @returns Named vector with the log likelihood for each of scars, dimples and detectedDimples that is not NULL
#'
#' @export
get.logLik.global.h <- function(h, field, scars = NULL, dimples = NULL, detectedDimples = NULL, l_x, l_y, tidspunkt,
                                localStat = "variance"){
  localStat <- match.arg(localStat, choices = c("variance", "meanSquare", "meanAbs"))
  pointSets <- collect.point.sets(scars = scars, dimples = dimples, detectedDimples = detectedDimples)
  if(any(lengths(pointSets) != length(tidspunkt))){
    stop("scars, dimples and detectedDimples must have one element per timestep in tidspunkt")
  }
  if(!is.numeric(tidspunkt)){
    stop("Experimental data, with tidspunkt as a list of (ensemble, timestep), is not supported yet")
  }

  Zglobal <- vapply(tidspunkt + h, function(t){
    switch(localStat, meanAbs = mean(abs(field[t, , ])), mean(field[t, , ]^2))
  }, numeric(1))

  return(vapply(pointSets, function(structures){
    fit.global.covar(Zglobal = Zglobal, structures = structures, l_x = l_x, l_y = l_y, returnParams = FALSE)
  }, numeric(1)))
}

#' Calculate the log likelihood for a specific radius and lag
#' @description
#' Compute the log likelihood(up to some constant) for a given local statistic, with a given radius and lag
#'
#' @param r Radius in pixels of the window used for the local statistic. Set to "infty" for the global model, where the
#' covariate is constant over each snapshot, which is calculated with get.logLik.global.h. If h is NULL, r is the radius
#' and lag encoded in one number, see make.input.matrix
#' @param h Lag, in timesteps, from the snapshots of the points to the snapshots of the field
#' @param field TxMxN array with the field the covariate is calculated from
#' @param scars List of scars at each timestep in tidspunkt, as returned by make.points.list, or NULL to leave them out
#' @param dimples List of dimples at each timestep in tidspunkt, as returned by make.points.list, or NULL to leave them
#' out
#' @param detectedDimples List of dimples detected from the surface elevation at each timestep in tidspunkt, as returned
#' by make.points.list, or NULL to leave them out
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param tidspunkt Vector with the timesteps of the points
#' @param localStat Statistic of the window used as covariate, "variance", "meanSquare" or "meanAbs". The model uses
#' log(Z), so the signed local mean can not be used. For the global model, "meanAbs" gives the mean absolute value of the
#' snapshot, and the other two give the mean square
#' @param radiiModulo Parameter for decoding r into radius and lag when h is NULL, see make.input.matrix
#' @param radiiForPrint Radii for which to print the radius, to follow the progress
#' @param startLag The radius is printed when the lag equals startLag
#'
#' @returns Named vector with the log likelihood for each of scars, dimples and detectedDimples that is not NULL
#'
#' @export
get.logLik.specific.r.h <- function(r, h = NULL, field, scars = NULL, dimples = NULL, detectedDimples = NULL, l_x, l_y,
                                    tidspunkt, localStat = "variance", radiiModulo = 1000L, radiiForPrint = -1,
                                    startLag = 0){
  localStat <- match.arg(localStat, choices = c("variance", "meanSquare", "meanAbs"))
  pointSets <- collect.point.sets(scars = scars, dimples = dimples, detectedDimples = detectedDimples)
  if(any(lengths(pointSets) != length(tidspunkt))){
    stop("scars, dimples and detectedDimples must have one element per timestep in tidspunkt")
  }
  if(is.null(h)){
    r_and_h_untransformed <- r
    r <- r_and_h_untransformed %% radiiModulo
    h <- (r_and_h_untransformed - r)/radiiModulo
    if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  }

  if(identical(r, "infty")){
    return(get.logLik.global.h(h = h, field = field, scars = scars, dimples = dimples, detectedDimples = detectedDimples,
                               l_x = l_x, l_y = l_y, tidspunkt = tidspunkt, localStat = localStat))
  }

  if(is.numeric(tidspunkt)){
    tidspunkt <- tidspunkt + h
  } else {
    for(t_ind in 1:length(tidspunkt)){
      tidspunkt[[t_ind]][2] <- tidspunkt[[t_ind]][2] + h
    }
  }

  if(!is.numeric(tidspunkt)){
    stop("Experimental data, with tidspunkt as a list of (ensemble, timestep), is not supported yet")
  }

  # drop = FALSE keeps the time dimension when there is only one snapshot
  field <- field[tidspunkt, , , drop = FALSE]
  Z <- array(data = NA_real_, dim = dim(field))
  if(!is.numeric(r) || r < 1){
    stop("r must be a radius of at least one pixel, or \"infty\" for the global model")
  }
  vindauga <- t(finn.vindauga.sirkel(radius = r))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
  calculateLocStat <- switch(localStat, variance = calculateLocVar_periodic,
                             meanSquare = calculateLocMeanSq_periodic, meanAbs = calculateLocMeanAbs_periodic)
  for(i in 1:length(tidspunkt)){
    # The C++ functions return the snapshot with the x index running fastest, which is how R fills Z[i, , ]
    Z[i, , ] <- calculateLocStat(inWin = vindauga, inShell = skall, insDiv = field[i, , ],
                                 inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2], dim(skall)[2]))
  }

  return(vapply(pointSets, function(structures){
    fit.given.covar(Z = Z, structures = structures, l_x = l_x, l_y = l_y, returnParams = FALSE)
  }, numeric(1)))
}


#' Make an array with the regression results for each choice of radii and lag
#' @description
#' Use this to calculate the log likelihood of the local model for scars, dimples and dimples detected from the surface
#' elevation, for every combination of radius and lag
#'
#' @param radiar Vector of radii in pixels. These must be whole numbers, since make.input.matrix stores them as integers
#' @param lags Vector of lags in timesteps
#' @param field TxMxN array with the field the covariate is calculated from
#' @param scars Complete list of scars for all timesteps, in the unit square, or NULL to leave them out
#' @param dimples Complete list of dimples for all timesteps, in the unit square, or NULL to leave them out
#' @param detectedDimples Complete list of dimples detected from the surface elevation for all timesteps, in the unit
#' square, or NULL to leave them out
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param tidspunkt Vector with the timesteps of the points
#' @param localStat Statistic of the window used as covariate, see get.logLik.specific.r.h
#' @param radiiModulo Parameter for encoding radius and lag in one number, see make.input.matrix
#' @param radiiForPrint Radii for which to print the radius, to follow the progress
#'
#' @returns List with one matrix for each of scars, dimples and detectedDimples that is not NULL, with one row per lag
#' and one column per radius
#'
#' @export
make.result.matrices <- function(radiar, lags, field, scars = NULL, dimples = NULL, detectedDimples = NULL, l_x, l_y,
                                 tidspunkt, localStat = "variance", radiiModulo = 1000L, radiiForPrint = -1){
  pointSets <- lapply(collect.point.sets(scars = scars, dimples = dimples, detectedDimples = detectedDimples),
                      make.points.list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)

  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  resArray <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.logLik.specific.r.h, h = NULL, field = field,
                    scars = pointSets$scars, dimples = pointSets$dimples, detectedDimples = pointSets$detectedDimples,
                    l_x = l_x, l_y = l_y, tidspunkt = tidspunkt, localStat = localStat, radiiModulo = radiiModulo,
                    radiiForPrint = radiiForPrint, startLag = lags[1])
  # apply leaves out the first dimension when there is only one point set
  resArray <- array(data = resArray, dim = c(length(pointSets), length(lags), length(radiar)))
  # matrix() keeps the result a matrix when there is only one lag or one radius
  res <- lapply(seq_along(pointSets), function(k){
    matrix(data = resArray[k, , ], nrow = length(lags), ncol = length(radiar), dimnames = list(lag = lags, radius = radiar))
  })
  return(stats::setNames(res, names(pointSets)))
}

#' Make vectors with the log likelihood of the global model for each lag
#' @description
#' Use this to calculate the log likelihood of the global model for scars, dimples and dimples detected from the
#' surface elevation, for every lag
#'
#' @param lags Vector of lags in timesteps
#' @param field TxMxN array with the field the covariate is calculated from
#' @param scars Complete list of scars for all timesteps, in the unit square, or NULL to leave them out
#' @param dimples Complete list of dimples for all timesteps, in the unit square, or NULL to leave them out
#' @param detectedDimples Complete list of dimples detected from the surface elevation for all timesteps, in the unit
#' square, or NULL to leave them out
#' @param l_x Window length in x direction
#' @param l_y Window length in y direction
#' @param tidspunkt Vector with the timesteps of the points
#' @param localStat Statistic of the snapshot used as covariate, see get.logLik.global.h
#'
#' @returns List with one vector for each of scars, dimples and detectedDimples that is not NULL, with one element per
#' lag
#'
#' @export
make.global.result.vectors <- function(lags, field, scars = NULL, dimples = NULL, detectedDimples = NULL, l_x, l_y,
                                       tidspunkt, localStat = "variance"){
  pointSets <- lapply(collect.point.sets(scars = scars, dimples = dimples, detectedDimples = detectedDimples),
                      make.points.list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)

  res <- vapply(X = lags, FUN = get.logLik.global.h, FUN.VALUE = numeric(length(pointSets)), field = field,
                scars = pointSets$scars, dimples = pointSets$dimples, detectedDimples = pointSets$detectedDimples,
                l_x = l_x, l_y = l_y, tidspunkt = tidspunkt, localStat = localStat)
  # vapply gives a vector instead of a matrix when there is only one point set
  res <- matrix(data = res, nrow = length(pointSets))
  return(stats::setNames(lapply(seq_along(pointSets), function(k){stats::setNames(res[k, ], lags)}), names(pointSets)))
}

#' Function for making a matrix containing information both on the lag and radius to be used for each index in the
#' result matrices
#' @description
#' Each element is radius + radiiModulo*lag, with one row per lag and one column per radius. This lets apply loop over
#' every combination of radius and lag, which get.logLik.specific.r.h and get.fit.specific.r.h decode when h is NULL.
#'
#' @param radiar Vector with all radii for local variance, whole numbers between 0 and radiiModulo
#' @param lags Vector with all lag values to use, whole numbers
#' @param radiiModulo parameter for how to encode lag and radius
#'
#' @export
make.input.matrix <- function(radiar, lags, radiiModulo = 1000L){
  if(any(radiar != round(radiar)) || any(radiar < 0) || max(radiar) >= radiiModulo){
    stop("radiar must be whole numbers between 0 and radiiModulo, otherwise the radii are not encoded correctly")
  }
  if(any(lags != round(lags))){
    stop("lags must be whole numbers, otherwise the lags are not encoded correctly")
  }
  res <- matrix(data = as.integer(radiar), ncol = length(radiar), nrow = length(lags), byrow = TRUE)
  for (h_i in 1:length(lags)) {
    res[h_i, ] <- res[h_i, ] + radiiModulo*lags[h_i]
  }
  return(res)
}

