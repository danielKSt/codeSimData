# library("rhdf5")
# library(codeSimData)
# library(spatstat)
#
# load(file = "/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/RE1000_WE10/simulatedPoints.RDa")
# load(file = "/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/RE1000_WE10/scales.RDa")
# field <- h5read(file = "/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/RE1000_WE10/RE1000_WE10_vorticity_z_0.000.mat", name = 'vortZ')
#
# # Test that basic components work: ----
#
# noFrames <- length(simulatedVortices)
# startTime <- 1L
# tidspunkt <- seq(from = startTime, to = noFrames, by = as.integer(ceiling(timescale)))
# rm(startTime, noFrames)
# r <- 5
# l <- dim(field)[2]/lengthscale
#
# scars_list <- make.points.list(points_list = simulatedScars, l = dim(field)[2]/lengthscale, tidspunkt = tidspunkt)
# dimples_list <- make.points.list(points_list = simulatedVortices, # l = 1, tidspunkt = tidspunkt)
#                                  l = dim(field)[2]/lengthscale, tidspunkt = tidspunkt)
# field_list <- lapply(X = tidspunkt, FUN = surf.transform.specific.t,
#                      field = field, l = l, vindauga = "infty")
#
# df <- spatstat.geom::hyperframe(scars = scars_list, dimples = dimples_list, field = field_list)
#
# mppm(scars ~ log(field), data = df)
#
#
# # Test of more complex functions work: ----
# fit_test <- get.fit.specific.r.h(scars_list = simulatedScars, dimples_list = simulatedVortices, l = dim(field)[2]/lengthscale,
#                      field = field, tidspunkt = seq(from = 1, to = length(simulatedVortices), by = ceiling(timescale)), r = 15, h = 0, resType = "fit")
#
# plot(fit_test$dimples.fit)
# summary(fit_test$dimples.fit)
#
# radiar <- seq(from = 1, to = 22, by = 4)
# lags <- seq(from = -50, to = 50, by = 10)
#
# inputMatrix <- make.input.matrix(radiar = seq(from = 1, to = 22, by = 4), lags = seq(from = -50, to = 50, by = 10))
#
# test_aic_mat <- make.aic.matrix(scars_list = simulatedScars,
#                 dimples_list = simulatedVortices,
#                 l = dim(field)[2]/lengthscale,
#                 field = field,
#                 tidspunkt = seq(from = 51, to = length(simulatedVortices)-51, by = ceiling(timescale)),
#                 radiar = seq(from = 1, to = 22, by = 4),
#                 lags = seq(from = -50, to = 50, by = 10),
#                 radiiForPrint = seq(from = 1, to = 22, by = 4))
#
# test_aic_fullWindow <- get.fit.specific.r.h(scars_list = simulatedScars, dimples_list = simulatedVortices, l = dim(field)[2]/lengthscale,
#                                             field = field, tidspunkt = seq(from = 1, to = length(simulatedVortices), by = ceiling(timescale)), r = "infty", h = 0, resType = "fit")

#' Function for preparing the points, given a list of all points
#' @param points_list List of all points in a unit square
#' @param l Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param tidspunkt Timesteps to use for list
#'
#' @export
make.points.list <- function(points_list, l, tidspunkt){
  res <- vector(mode = "list", length = length(tidspunkt))
  for(i in 1:length(tidspunkt)){
    t <- tidspunkt[i]
    if(is.numeric(points_list[[t]])){
      res[[i]] <- spatstat.geom::ppp(x = c(), y = c(), c(0,l), c(0,l))
    } else {
      res[[i]] <- spatstat.geom::ppp(x = points_list[[t]]$x*l, y = points_list[[t]]$y*l, c(0,l), c(0,l))
    }
  }
  return(res)
}

#' Function for preparing the points, given a list of all points
#' @param field TxMxN array with field to be transformed into covariate
#' @param l Side length of observation window
#' @param t Timestep to use of array
#' @param vindauga window to use for the calculation of local variance
#' @param skall shell of window to be used for calculation of local variance
#'
#' @export
surf.transform.specific.t <- function(t, field, l, vindauga = NULL, skall = NULL){
  if(is.null(vindauga)){
    Z <- spatstat.geom::im(t(field[t, , ]),
                           seq(from = 0, to = l, length = dim(field)[2]),
                           seq(from = 0, to = l, length = dim(field)[3]))
    return(Z)
  } else if(is.character(vindauga)) {
    ms_field <- mean(field[t, , ]^2)
    return(ms_field)
  } else {
    locVar <- matrix(data = calculateLocVar_periodic(inWin = vindauga, inShell = skall,
                                       insDiv = field[t, , ],
                                       inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2], dim(skall)[2])),
                     nrow = dim(field)[2], ncol = dim(field)[3])

    Z <- spatstat.geom::im(t(locVar),
                           seq(from = 0, to = l, length = dim(field)[2]),
                           seq(from = 0, to = l, length = dim(field)[3]))
    return(Z)
  }
}


#' Function for preparing the points, given a list of all points
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l Side length of window
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param r radius
#' @param h lag
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#' @param startLag Integer, if r is in radiiForPrint, and h == startLag, we print r
#' @param resType What type of output is desired? Set to "aic" to get the aic of the fitted model. Set to "fit" for the model object. Set to "logLik" for the log-likelihood
#'
#' @export
get.fit.specific.r.h <- function(scars_list, dimples_list, l, field, tidspunkt, r, h = 0, radiiModulo = 1000L, radiiForPrint = -1, startLag = 0, resType = "logLik"){
  if(is.null(h)){
    r_and_h_untransformed <- r
    r <- r_and_h_untransformed %% radiiModulo
    h <- (r_and_h_untransformed - r)/radiiModulo
    if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  }
  scars_list <- make.points.list(points_list = scars_list, l = l, tidspunkt = tidspunkt)
  dimples_list <- make.points.list(points_list = dimples_list, l = l, tidspunkt = tidspunkt)
  tidspunkt <- tidspunkt + h

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
                       field = field, l = l, vindauga = vindauga, skall = skall)

  df <- spatstat.geom::hyperframe(scars = scars_list, dimples = dimples_list, field = field_list)
  scars_log.fit <- spatstat.model::mppm(formula = scars ~ log(field), data = df)
  scars.fit <- spatstat.model::mppm(formula = scars ~ field, data = df)
  dimples_log.fit <- spatstat.model::mppm(formula = dimples ~ log(field), data = df)
  dimples.fit <- spatstat.model::mppm(formula = dimples ~ field, data = df)
  if(resType == "logLik"){
    return(c(stats::logLik(scars_log.fit), stats::logLik(dimples_log.fit), stats::logLik(scars.fit), stats::logLik(dimples.fit)))
  } else if(resType == "aic"){
    return(c(stats::AIC(scars_log.fit), stats::AIC(dimples_log.fit), stats::AIC(scars.fit), stats::AIC(dimples.fit)))
  } else if (resType == "fit"){
    return(list(scars.fit = scars_log.fit, dimples.fit = dimples_log.fit, scars_noLog.fit = scars.fit, dimples_noLog.fit = dimples.fit))
  }
}

#' Function to make an aic matrix using spatstats mppm function
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l Side length of window
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param radiar radiar
#' @param lags lags
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#'
#' @export
make.aic.matrix <- function(scars_list, dimples_list, l, field, tidspunkt, radiar, lags, radiiModulo = 1000L, radiiForPrint){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.fit.specific.r.h, scars_list = scars_list, l = l,
                       dimples_list = dimples_list, field = field, tidspunkt = tidspunkt, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, h = NULL, resType = "aic")

  res.scars <- aicMatrices[1, , ]
  res.dimples <- aicMatrices[2, , ]
  res.scars_logless <- aicMatrices[3, , ]
  res.dimples_logless <- aicMatrices[4, , ]

  return(list(scars = res.scars, dimples = res.dimples, scars.noLog = res.scars_logless, dimples.noLog = res.dimples_logless))
}

#' Function to make an aic matrix using spatstats mppm function
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l Side length of window
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param radiar radiar
#' @param lags lags
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#'
#' @export
make.logLik.matrix <- function(scars_list, dimples_list, l, field, tidspunkt, radiar, lags, radiiModulo = 1000L, radiiForPrint){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  logLikMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.fit.specific.r.h, scars_list = scars_list, l = l,
                       dimples_list = dimples_list, field = field, tidspunkt = tidspunkt, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, h = NULL, resType = "logLik")

  res.scars <- logLikMatrices[1, , ]
  res.dimples <- logLikMatrices[2, , ]
  res.scars_logless <- logLikMatrices[3, , ]
  res.dimples_logless <- logLikMatrices[4, , ]

  return(list(scars = res.scars, dimples = res.dimples, scars.noLog = res.scars_logless, dimples.noLog = res.dimples_logless))
}

