# Fitting the regression model with mppm from spatstat. The exact likelihood in result_matrices.R does not need
# spatstat, but this is still the only fitting that supports experimental data.

#' Function for preparing the points, given a list of all points
#' @param points_list List of all points in a unit square
#' @param l_x Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param l_y Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param tidspunkt Timesteps to use for list
#'
#' @export
make.points.list.spatstat <- function(points_list, l_x, l_y, tidspunkt){
  res <- vector(mode = "list", length = length(tidspunkt))
  if(is.numeric(tidspunkt)){
    for(i in 1:length(tidspunkt)){
      t <- tidspunkt[i]
      if(is.empty.snapshot(points_list[[t]])){
        res[[i]] <- spatstat.geom::ppp(x = c(), y = c(), c(0,l_x), c(0,l_y))
      } else {
        res[[i]] <- spatstat.geom::ppp(x = points_list[[t]]$x*l_x, y = points_list[[t]]$y*l_y, c(0,l_x), c(0,l_y))
      }
    }
  } else {
    for(i in 1:length(tidspunkt)){
      ensemble <- tidspunkt[[i]][1]
      t <- tidspunkt[[i]][2]
      if(is.empty.snapshot(points_list[[ensemble]][[t]])){
        res[[i]] <- spatstat.geom::ppp(x = c(), y = c(), c(0,l_x), c(0,l_y))
      } else {
        res[[i]] <- spatstat.geom::ppp(x = points_list[[ensemble]][[t]]$x*l_x, y = points_list[[ensemble]][[t]]$y*l_y, c(0,l_x), c(0,l_y))
      }
    }
  }
  return(res)
}

#' Transform the scalar field with a local statistic (by default the local variance) for a given window, or the full window MSE
#' @param field TxMxN array with field to be transformed into covariate
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param t Timestep to use of array, or a vector of the form
#' @param vindauga window to use for the calculation of the local statistic
#' @param skall shell of window to be used for calculation of the local statistic
#' @param localStat Statistic to calculate over the window. "variance" for the local variance, "meanSquare" for the sum
#' of the squared values divided by the number of pixels in the window, "meanAbs" for the sum of the absolute values
#' divided by the number of pixels in the window, and "mean" for the sum of the values divided by the number of pixels
#' in the window. If vindauga is a character, "mean" and "meanAbs" give the mean and the mean absolute value of the full
#' field, and the other two give the full window MSE
#'
#' @export
surf.transform.specific.t <- function(t, field, l_x, l_y, vindauga = NULL, skall = NULL, localStat = "variance"){
  localStat <- match.arg(localStat, choices = c("variance", "meanSquare", "meanAbs", "mean"))
  periodic <- TRUE
  if(length(t) == 2){
    field <- field[t[1], , , ]
    t <- t[2]
    periodic <- FALSE
  }
  # Pixel k covers ((k-1)*l/N, k*l/N), matching the point coordinates from explore.vortex and explore.scar
  if(is.null(vindauga)){
    Z <- spatstat.geom::im(t(field[t, , ]), xrange = c(0, l_x), yrange = c(0, l_y))
    return(Z)
  } else if(is.character(vindauga)) {
    if(localStat == "mean"){
      ms_field <- mean(field[t, , ])
    } else if(localStat == "meanAbs"){
      ms_field <- mean(abs(field[t, , ]))
    } else {
      ms_field <- mean(field[t, , ]^2)
    }
    Z <- spatstat.geom::im(matrix(data = ms_field, nrow = dim(field)[3], ncol = dim(field)[2]),
                           xrange = c(0, l_x), yrange = c(0, l_y))
    return(Z)
  } else {
    if(!periodic){
      calculateLocStat <- switch(localStat, variance = calculateLocVar_physical,
                                 meanSquare = calculateLocMeanSq_physical, meanAbs = calculateLocMeanAbs_physical,
                                 mean = calculateLocMean_physical)
    } else {
      calculateLocStat <- switch(localStat, variance = calculateLocVar_periodic,
                                 meanSquare = calculateLocMeanSq_periodic, meanAbs = calculateLocMeanAbs_periodic,
                                 mean = calculateLocMean_periodic)
    }
    locStat <- matrix(data = calculateLocStat(inWin = vindauga, inShell = skall,
                                              insDiv = field[t, , ],
                                              inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2], dim(skall)[2])),
                      nrow = dim(field)[2], ncol = dim(field)[3])

    Z <- spatstat.geom::im(t(locStat), xrange = c(0, l_x), yrange = c(0, l_y))
    return(Z)
  }
}


#' Fit the model with mppm for a specific radius and lag
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param r radius. Set to 0 to use the field itself as covariate, which then has to be positive, or to "infty" for the
#' global model
#' @param h lag
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#' @param startLag Integer, if r is in radiiForPrint, and h == startLag, we print r
#' @param resType What type of output is desired? Set to "aic" to get the aic of the fitted model. Set to "fit" for the model object. Set to "logLik" for the log-likelihood
#' @param printProgressComplete Set to TRUE if you want detailed updates on the progress
#' @param fitLogless Set to TRUE if you want to fit the model directly to the transformed field as well as the logarithm
#' @param localStat Statistic of the window used as covariate, "variance", "meanSquare", "meanAbs" or "mean", see surf.transform.specific.t.
#' The log fits need a positive covariate, so "mean" only works for positive fields
#' @param pixelQuadrature Set to TRUE to use one dummy point per pixel centre. The covariate is constant within each
#' pixel, so this makes the likelihood exact. Set to FALSE to use the default quadrature of mppm, which is faster, but
#' only approximates the integral of the intensity
#'
#' @returns For resType "logLik" and "aic", a named vector with elements scars and dimples, and also scars.noLog and
#' dimples.noLog if fitLogless is TRUE. For resType "fit", a list with the fitted models
#'
#' @export
get.fit.specific.r.h <- function(scars_list, dimples_list, l_x, l_y, field, tidspunkt, r, h = 0, radiiModulo = 1000L,
                                 radiiForPrint = -1, startLag = 0,
                                 resType = "logLik", printProgressComplete = FALSE, fitLogless = FALSE,
                                 localStat = "variance", pixelQuadrature = TRUE){
  resType <- match.arg(resType, choices = c("logLik", "aic", "fit"))
  if(is.null(h)){
    r_and_h_untransformed <- r
    r <- r_and_h_untransformed %% radiiModulo
    h <- (r_and_h_untransformed - r)/radiiModulo
    if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  }
  if(printProgressComplete){
    time_prev <- Sys.time()
  }
  scars_list <- make.points.list.spatstat(points_list = scars_list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)
  dimples_list <- make.points.list.spatstat(points_list = dimples_list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)
  if(printProgressComplete){
    time_curr <- Sys.time()
    print(paste("Time on preparing scars and dimples for r = ", r, ", h = ", h, ":  ",
                difftime(time_curr-time_prev, units = "mins"), "mins" , sep = ""))
    time_prev <- time_curr
  }
  if(is.numeric(tidspunkt)){
    tidspunkt <- tidspunkt + h
  } else {
    for(t_ind in 1:length(tidspunkt)){
      tidspunkt[[t_ind]][2] <- tidspunkt[[t_ind]][2] + h
    }
  }
  if(printProgressComplete){
    time_curr <- Sys.time()
    print(paste("Time on applying lag for r = ", r, ", h = ", h, ":  ",
                difftime(time_curr-time_prev, units = "mins"), "mins" , sep = ""))
    time_prev <- time_curr
  }

  if(r == 0){
    vindauga <- NULL
    skall <- NULL
  } else if(is.numeric(r)) {
    vindauga <- t(finn.vindauga.sirkel(radius = r))
    skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
  } else {
    vindauga <- "infty"
    skall <- NULL
  }
  field_list <- lapply(X = tidspunkt, FUN = surf.transform.specific.t,
                       field = field, l_x = l_x, l_y = l_y, vindauga = vindauga, skall = skall, localStat = localStat)
  if(printProgressComplete){
    time_curr <- Sys.time()
    print(paste("Time on transforming spatial field for r = ", r, ", h = ", h, ":  ",
                difftime(time_curr-time_prev, units = "mins"), "mins" , sep = ""))
    time_prev <- time_curr
  }
  # The field itself (r = 0) and the local mean can be negative, and log(field) is then undefined
  if(any(vapply(field_list, function(Z){min(Z$v, na.rm = TRUE) <= 0}, logical(1)))){
    stop("The covariate must be positive, since the model uses log(field). For a signed field, use r >= 1 with ",
         "localStat \"variance\", \"meanSquare\" or \"meanAbs\"")
  }
  if(pixelQuadrature){
    to.pixel.quad <- function(X, Z){spatstat.geom::pixelquad(X, W = spatstat.geom::as.owin(Z))}
    scars_list <- mapply(to.pixel.quad, scars_list, field_list, SIMPLIFY = FALSE)
    dimples_list <- mapply(to.pixel.quad, dimples_list, field_list, SIMPLIFY = FALSE)
  }

  df <- spatstat.geom::hyperframe(scars = scars_list, dimples = dimples_list, field = field_list)
  scars_log.fit <- spatstat.model::mppm(formula = scars ~ log(field), data = df)
  dimples_log.fit <- spatstat.model::mppm(formula = dimples ~ log(field), data = df)
  if(fitLogless){
    dimples.fit <- spatstat.model::mppm(formula = dimples ~ field, data = df)
    scars.fit <- spatstat.model::mppm(formula = scars ~ field, data = df)
  }
  if(printProgressComplete){
    time_curr <- Sys.time()
    print(paste("Time on fitting model for r = ", r, ", h = ", h, ":  ",
                difftime(time_curr-time_prev, units = "mins"), "mins", sep = ""))
    time_prev <- time_curr
  }
  if (resType == "fit"){
    if(fitLogless){
      return(list(scars.fit = scars_log.fit, dimples.fit = dimples_log.fit, scars_noLog.fit = scars.fit, dimples_noLog.fit = dimples.fit))
    } else {
      return(list(scars.fit = scars_log.fit, dimples.fit = dimples_log.fit))
    }
  }
  measure <- switch(resType, logLik = function(fit){as.numeric(stats::logLik(fit))}, aic = stats::AIC)
  res <- c(scars = measure(scars_log.fit), dimples = measure(dimples_log.fit))
  if(fitLogless){
    res <- c(res, scars.noLog = measure(scars.fit), dimples.noLog = measure(dimples.fit))
  }
  return(res)
}

#' Function to make an aic matrix using spatstats mppm function
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param radiar radiar
#' @param lags lags
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#' @param printProgressComplete Set to TRUE if you want detailed updates on the progress
#' @param fitLogless Set to TRUE if you want to fit the model directly to the transformed field as well as the logarithm
#' @param localStat Statistic of the window used as covariate, "variance", "meanSquare", "meanAbs" or "mean", see surf.transform.specific.t.
#' The log fits need a positive covariate, so "mean" only works for positive fields
#' @param pixelQuadrature Set to TRUE to use one dummy point per pixel centre, which makes the likelihood exact, see
#' get.fit.specific.r.h
#'
#' @returns List with the matrices scars and dimples, and also scars.noLog and dimples.noLog if fitLogless is TRUE.
#' Each matrix has one row per lag and one column per radius
#'
#' @export
make.logLik.matrix <- function(scars_list, dimples_list, l_x, l_y, field, tidspunkt, radiar, lags, radiiModulo = 1000L,
                               radiiForPrint, printProgressComplete = FALSE, fitLogless = FALSE, localStat = "variance",
                               pixelQuadrature = TRUE){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  logLikMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.fit.specific.r.h, scars_list = scars_list, l_x = l_x, l_y = l_y,
                       dimples_list = dimples_list, field = field, tidspunkt = tidspunkt, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, h = NULL, resType = "logLik",
                       printProgressComplete = printProgressComplete, fitLogless = fitLogless, localStat = localStat,
                       pixelQuadrature = pixelQuadrature)

  # matrix() keeps the result a matrix when there is only one lag or one radius
  as.result.matrix <- function(name){
    matrix(data = logLikMatrices[name, , ], nrow = length(lags), ncol = length(radiar),
           dimnames = list(lag = lags, radius = radiar))
  }
  res <- list(scars = as.result.matrix("scars"), dimples = as.result.matrix("dimples"))
  if(fitLogless){
    res$scars.noLog <- as.result.matrix("scars.noLog")
    res$dimples.noLog <- as.result.matrix("dimples.noLog")
  }
  return(res)
}
