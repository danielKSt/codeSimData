# library("rhdf5")
# library(codeSimData)
# library(spatstat)
#
# load("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/structures/pointsP50.RDa")
# load("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/scales/p50.RDa")
# hDiv <- h5read(file = "/Volumes/LaCie/Daniel/Experimental/p50_h390.h5", name = "hDiv")
#
# field <- hDiv[1, , , ]
# image(field[1, , ])
# image(hDiv[1, 1, , ])
# tidspunkt <- list(c(1, 23), c(1, 50), c(2, 10), c(2, 50))
#
# a <- make.points.list(points_list = dimples, l_x = dim(hDiv)[3]/lengthscale, l_y = dim(hDiv)[4]/lengthscale, tidspunkt = tidspunkt)
# plot(a[[1]])
#
# b <- surf.transform.specific.t(t = c(1, 1), field = hDiv, l_x = dim(hDiv)[3]/lengthscale, l_y = dim(hDiv)[4]/lengthscale, vindauga = NULL)
# b <- surf.transform.specific.t(t = c(1, 1), field = hDiv, l_x = dim(hDiv)[3]/lengthscale, l_y = dim(hDiv)[4]/lengthscale, vindauga = "NULL")
#
# vindauga <- t(finn.vindauga.sirkel(radius = 5))
# skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
# b <- surf.transform.specific.t(t = c(1, 1), field = hDiv, l_x = dim(hDiv)[3]/lengthscale, l_y = dim(hDiv)[4]/lengthscale, vindauga = vindauga, skall = skall)
# plot(b)

#' Function for preparing the points, given a list of all points
#' @param points_list List of all points in a unit square
#' @param l_x Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param l_x Sidelength of box, both for scaling points up to correct size, and to create owin objects
#' @param tidspunkt Timesteps to use for list
#'
#' @export
make.points.list <- function(points_list, l_x, l_y, tidspunkt){
  res <- vector(mode = "list", length = length(tidspunkt))
  if(is.numeric(tidspunkt)){
    for(i in 1:length(tidspunkt)){
      t <- tidspunkt[i]
      if(is.numeric(points_list[[t]])){
        res[[i]] <- spatstat.geom::ppp(x = c(), y = c(), c(0,l_x), c(0,l_y))
      } else {
        res[[i]] <- spatstat.geom::ppp(x = points_list[[t]]$x*l_x, y = points_list[[t]]$y*l_y, c(0,l_x), c(0,l_y))
      }
    }
  } else {
    for(i in 1:length(tidspunkt)){
      ensemble <- tidspunkt[[i]][1]
      t <- tidspunkt[[i]][2]
      if(is.numeric(points_list[[ensemble]][[t]])){
        res[[i]] <- spatstat.geom::ppp(x = c(), y = c(), c(0,l_x), c(0,l_y))
      } else {
        res[[i]] <- spatstat.geom::ppp(x = points_list[[ensemble]][[t]]$x*l_x, y = points_list[[ensemble]][[t]]$y*l_y, c(0,l_x), c(0,l_y))
      }
    }
  }
  return(res)
}

#' Transform the scalar field with either the local variance for a given window, or the full window MSE
#' @param field TxMxN array with field to be transformed into covariate
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param t Timestep to use of array, or a vector of the form
#' @param vindauga window to use for the calculation of local variance
#' @param skall shell of window to be used for calculation of local variance
#'
#' @export
surf.transform.specific.t <- function(t, field, l_x, l_y, vindauga = NULL, skall = NULL){
  periodic <- TRUE
  if(length(t) == 2){
    field <- field[t[1], , , ]
    t <- t[2]
    periodic <- FALSE
  }
  if(is.null(vindauga)){
    Z <- spatstat.geom::im(t(field[t, , ]),
                           seq(from = 0, to = l_x, length = dim(field)[2]),
                           seq(from = 0, to = l_y, length = dim(field)[3]))
    return(Z)
  } else if(is.character(vindauga)) {
    ms_field <- mean(field[t, , ]^2)
    return(ms_field)
  } else {
    if(!periodic){
      locVar <- matrix(data = calculateLocVar_physical(inWin = vindauga, inShell = skall,
                                                       insDiv = field[t, , ],
                                                       inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2], dim(skall)[2])),
                       nrow = dim(field)[2], ncol = dim(field)[3])

      Z <- spatstat.geom::im(t(locVar),
                             seq(from = 0, to = l_x, length = dim(field)[2]),
                             seq(from = 0, to = l_y, length = dim(field)[3]))
      return(Z)
    } else {
      locVar <- matrix(data = calculateLocVar_periodic(inWin = vindauga, inShell = skall,
                                                       insDiv = field[t, , ],
                                                       inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2], dim(skall)[2])),
                       nrow = dim(field)[2], ncol = dim(field)[3])

      Z <- spatstat.geom::im(t(locVar),
                             seq(from = 0, to = l_x, length = dim(field)[2]),
                             seq(from = 0, to = l_y, length = dim(field)[3]))
      return(Z)
    }
  }
}


#' Function for preparing the points, given a list of all points
#' @param scars_list Complete list of scars for all timesteps
#' @param dimples_list Complete list of dimples for all timesteps
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
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
get.fit.specific.r.h <- function(scars_list, dimples_list, l_x, l_y, field, tidspunkt, r, h = 0, radiiModulo = 1000L, radiiForPrint = -1, startLag = 0, resType = "logLik"){
  if(is.null(h)){
    r_and_h_untransformed <- r
    r <- r_and_h_untransformed %% radiiModulo
    h <- (r_and_h_untransformed - r)/radiiModulo
    if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  }
  scars_list <- make.points.list(points_list = scars_list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)
  dimples_list <- make.points.list(points_list = dimples_list, l_x = l_x, l_y = l_y, tidspunkt = tidspunkt)
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
                       field = field, l_x = l_x, l_y = l_y, vindauga = vindauga, skall = skall)

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
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param radiar radiar
#' @param lags lags
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#'
#' @export
make.aic.matrix <- function(scars_list, dimples_list, l_x, l_y, field, tidspunkt, radiar, lags, radiiModulo = 1000L, radiiForPrint){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.fit.specific.r.h, scars_list = scars_list, l_x = l_x, l_y = l_y,
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
#' @param l_x Side length of observation window in x direction
#' @param l_y Side length of observation window in y direction
#' @param field TxMxN array with surface field
#' @param tidspunkt vector with timesteps to use for inference
#' @param radiar radiar
#' @param lags lags
#' @param radiiModulo Parameter for the case where we use a single number for both r and h
#' @param radiiForPrint Vector of all radii where we should print the radius
#'
#' @export
make.logLik.matrix <- function(scars_list, dimples_list, l_x, l_y, field, tidspunkt, radiar, lags, radiiModulo = 1000L, radiiForPrint){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  logLikMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = get.fit.specific.r.h, scars_list = scars_list, l_x = l_x, l_y = l_y,
                       dimples_list = dimples_list, field = field, tidspunkt = tidspunkt, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, h = NULL, resType = "logLik")

  res.scars <- logLikMatrices[1, , ]
  res.dimples <- logLikMatrices[2, , ]
  res.scars_logless <- logLikMatrices[3, , ]
  res.dimples_logless <- logLikMatrices[4, , ]

  return(list(scars = res.scars, dimples = res.dimples, scars.noLog = res.scars_logless, dimples.noLog = res.dimples_logless))
}

