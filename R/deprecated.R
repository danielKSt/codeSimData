# Deprecated: analysis of DNS data by stacking the snapshots into one tall image and fitting with ppm from spatstat,
# with and without a shrunk observation window. Use make.result.matrices or make.logLik.matrix instead.

# Records whether the deprecation warning has been given in this session
deprecationState <- new.env(parent = emptyenv())
deprecationState$warned <- FALSE

#' Warn that the stacked-image analysis of DNS data is deprecated
#' @description
#' The warning is only given the first time a deprecated function is called in a session, since the functions are
#' called once for every radius and lag. Set the option codeSimData.silenceDeprecation to TRUE to turn it off.
#'
#' @param fun Name of the deprecated function that was called
warn.stacked.deprecated <- function(fun){
  if(deprecationState$warned || isTRUE(getOption("codeSimData.silenceDeprecation"))){
    return(invisible(NULL))
  }
  deprecationState$warned <- TRUE
  warning(fun, " is deprecated, along with the rest of the stacked-image analysis of DNS data. ",
          "Use make.result.matrices or make.logLik.matrix instead. This warning is only shown once per session.",
          call. = FALSE)
}

#' Deprecated: function for preparing the points, given a list of all points
#' @param indices vector of indices/times to use points from
#' @param points The complete list of points at all times
#'
#' @keywords internal
#' @export
prepare.points.dns <- function(indices, points){
  warn.stacked.deprecated("prepare.points.dns")
  t <- 1
  firstPointNotFound <- TRUE
  while(firstPointNotFound){
    if (!is.empty.snapshot(points[[indices[t]]])) {
      temp <- rbind(points[[indices[t]]])

      for(p in 1:length(temp$x)){
        if(temp$x[p] < 0){
          temp$x[p] <- temp$x[p] + 1
        }
        if(temp$x[p] > 1){
          temp$x[p] <- temp$x[p] - 1
        }
        if(temp$y[p] < 0){
          temp$y[p] <- temp$y[p] + 1
        }
        if(temp$y[p] > 1){
          temp$y[p] <- temp$y[p] - 1
        }
      }

      if(dim(temp)[1] > 0){
        firstPointNotFound <- FALSE
        l <- dim(temp)[1]
        pointsForAnalysis <- rbind(temp+c(rep(0, l), rep(t-1, l)))
        break
      }
    }
    t <- t + 1
    if(t > length(indices)){
      print("No points at given times.")
      break
    }
  }

  timeOfFirstPoint <- t
  for (t in (timeOfFirstPoint+1):length(indices)) {
    if(!is.empty.snapshot(points[[indices[t]]])){
      temp <- rbind(points[[indices[t]]])

      for(p in 1:length(temp$x)){
        if(temp$x[p] < 0){
          temp$x[p] <- temp$x[p] + 1
        }
        if(temp$x[p] > 1){
          temp$x[p] <- temp$x[p] - 1
        }
        if(temp$y[p] < 0){
          temp$y[p] <- temp$y[p] + 1
        }
        if(temp$y[p] > 1){
          temp$y[p] <- temp$y[p] - 1
        }
      }

      if(dim(temp)[1] > 0){
        l <- dim(temp)[1]
        pointsForAnalysis <- rbind(pointsForAnalysis, temp+c(rep(0, l), rep(t-1, l)))
      }
    }
  }


  points.pp <- spatstat.geom::ppp(x = pointsForAnalysis$x, y = pointsForAnalysis$y, c(0,1), c(0,length(indices)))
  rm(pointsForAnalysis, temp)

  return(points.pp)
}

#' Deprecated: function for creating local variance covariates for the desired indices for DNS data
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param r radius to use in calculating local variance
#' @param h lag
#'
#' @keywords internal
#' @export
combine.covariates.dns <- function(hDiv, indices, r, h){
  warn.stacked.deprecated("combine.covariates.dns")
  numberOfPictures <- length(indices)

  xdim <- dim(hDiv)[2]
  ydim <- dim(hDiv)[3]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  vindauga <- t(finn.vindauga.sirkel(radius = r))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for (t in 1:numberOfPictures) {
    insDiv <- hDiv[indices[t]+h, , ]
    covariatesToCombine[[t]] <- calculateLocVar_periodic(inWin = vindauga, inShell = skall,
                                                         insDiv = insDiv,
                                                         inDims = c(xdim, ydim, dim(vindauga)[2], dim(skall)[2]))
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim, ydim = ydim)

  Z <- spatstat.geom::im(t(combinedCovariates), xrange = c(0, 1), yrange = c(0, numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, hDiv, r, h)

  return(Z)
}

#' Deprecated: function for creating beta squared covariates for the desired indices for DNS data
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param h lag
#'
#' @keywords internal
#' @export
combine.covariates.betasq.dns <- function(hDiv, indices, h){
  warn.stacked.deprecated("combine.covariates.betasq.dns")
  numberOfPictures <- length(indices)

  xdim <- dim(hDiv)[2]
  ydim <- dim(hDiv)[3]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (t in 1:numberOfPictures) {
    insDiv <- hDiv[indices[t]+h, , ]
    covariatesToCombine[[t]] <- sum(insDiv^2) * matrix(data = 1, nrow = xdim, ncol = ydim)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim, ydim = ydim)

  Z <- spatstat.geom::im(t(combinedCovariates), xrange = c(0, 1), yrange = c(0, numberOfPictures))
  rm(covariatesToCombine, xdim, ydim, numberOfPictures, hDiv, h)

  return(Z)
}

#' Deprecated: function for making a new AIC matrix given DNS data for a given radius and lag
#' @param r_and_h_untransformed r_and_h_combined in single number
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param hDiv horisontal divergence
#' @param indices which indices to use
#' @param radiiForPrint which radii to print
#' @param startLag First lag, used to decide when to print the radius
#' @param radiiModulo Parameter for how to decode the r_and_h_untransformed to r and h
#' @param transformation Should the covariate be transformed, if value is "log" then a log transformation is applied
#'
#' @keywords internal
#' @export
aic.given.radius.and.lag <- function(r_and_h_untransformed, vortices.pp, scars.pp, hDiv,
                                    indices, radiiForPrint, startLag, radiiModulo = 1000L, transformation = "identity"){
  warn.stacked.deprecated("aic.given.radius.and.lag")
  r <- r_and_h_untransformed %% radiiModulo
  h <- (r_and_h_untransformed - r)/radiiModulo
  if((r %in% radiiForPrint) && (h == startLag)) {print(r)}

  Z <- combine.covariates.dns(hDiv = hDiv, indices = indices, r = r, h = h)
  if(transformation == "log"){Z <- log(Z)}
  vort.fit <- spatstat.model::ppm(vortices.pp ~ Z)
  scar.fit <- spatstat.model::ppm(scars.pp ~ Z)
  if(transformation == "both"){
    vort.fit.log <- spatstat.model::ppm(vortices.pp ~ log(Z))
    scar.fit.log <- spatstat.model::ppm(scars.pp ~ log(Z))
    rm(Z, r, h)
    return(c(stats::AIC(vort.fit), stats::AIC(scar.fit), stats::AIC(vort.fit.log), stats::AIC(scar.fit.log)))
  }

  rm(Z, r, h)
  return(c(stats::AIC(vort.fit), stats::AIC(scar.fit)))
}

#' Function for scaling points from a reduced observation window into a "full" observation window
#' @param points dataframe of points
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
scale_points_shrunk <- function(points, xlim, ylim){
  points$x <- (points$x-xlim[1])/(xlim[2]-xlim[1])
  points$y <- (points$y-ylim[1])/(ylim[2]-ylim[1])
  return(points)
}


#' Function for checking if any of the points are in the shrunk observation window
#' @param points dataframe of points
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
point.in.shrunk <- function(points, xlim, ylim){
  for(i in 1:nrow(points)){
    inShrunk <- ((xlim[1] <= points$x[i]) && (xlim[2] >= points$x[i]) && (ylim[1] <= points$y[i]) && (ylim[2] >= points$y[i]))
    if(inShrunk){
      break
    }
  }
  return(inShrunk)
}


#' Function for checking if any of the points are in the shrunk observation window
#' @param points dataframe of points
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
find.points.in.shrunk <- function(points, xlim, ylim){
  inShrunk <- c(1:nrow(points))
  for(i in 1:nrow(points)){
    inShrunk[i] <- ((xlim[1] <= points$x[i]) && (xlim[2] >= points$x[i]) && (ylim[1] <= points$y[i]) && (ylim[2] >= points$y[i]))
  }
  return(which(inShrunk == 1))
}

#' Deprecated: function for preparing the points, given a list of all points
#' @param indices vector of indices/times to use points from
#' @param points The complete list of points at all times
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @keywords internal
#' @export
prepare.points.dns.shrunk <- function(indices, points, xlim, ylim){
  warn.stacked.deprecated("prepare.points.dns.shrunk")
  t <- 1
  firstPointNotFound <- TRUE
  while(firstPointNotFound){
    if (!is.empty.snapshot(points[[indices[t]]])) {
      if (point.in.shrunk(points = points[[indices[t]]], xlim = xlim, ylim = ylim)){
        temp <- points[[indices[t]]]
        temp <- temp[find.points.in.shrunk(points = temp, xlim = xlim, ylim = ylim), ]
        temp <- scale_points_shrunk(points = temp, xlim = xlim, ylim = ylim)
        temp <- rbind(temp)
        if(dim(temp)[1] > 0){
          firstPointNotFound <- FALSE
          l <- dim(temp)[1]
          pointsForAnalysis <- rbind(temp+c(rep(0, l), rep(t-1, l)))
          break
        }
      }
    }
    t <- t + 1
    if(t > length(indices)){
      print("No points at given times.")
      break
    }
  }
  timeOfFirstPoint <- t
  for (t in (timeOfFirstPoint+1):length(indices)) {
    if(!is.empty.snapshot(points[[indices[t]]])){
      if(point.in.shrunk(points = points[[indices[t]]], xlim = xlim, ylim = ylim)){
        temp <- points[[indices[t]]]
        temp <- temp[find.points.in.shrunk(points = temp, xlim = xlim, ylim = ylim), ]
        temp <- scale_points_shrunk(points = temp, xlim = xlim, ylim = ylim)
        temp <- rbind(temp)
        if(dim(temp)[1] > 0){
          l <- dim(temp)[1]
          pointsForAnalysis <- rbind(pointsForAnalysis, temp+c(rep(0, l), rep(t-1, l)))
        }
      }
    }
  }
  points.pp <- spatstat.geom::ppp(x = pointsForAnalysis$x, y = pointsForAnalysis$y, c(0,1), c(0,length(indices)))
  rm(pointsForAnalysis, temp)
  return(points.pp)
}

#' Deprecated: function for creating local variance covariates for the desired indices for DNS data
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param r radius to use in calculating local variance
#' @param h lag
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @keywords internal
#' @export
combine.covariates.dns.shrunk <- function(sDiv, indices, r, h, xlim, ylim){
  warn.stacked.deprecated("combine.covariates.dns.shrunk")
  numberOfPictures <- length(indices)

  xdim <- dim(sDiv)[2]
  xdimShrunk <- xlim[2] - xlim[1] + 1
  ydim <- dim(sDiv)[3]
  ydimShrunk <- ylim[2] - ylim[1] + 1

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  vindauga <- t(finn.vindauga.sirkel(radius = r))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for (t in 1:numberOfPictures) {
    insDiv <- sDiv[indices[t]+h, , ]
    temp <- calculateLocVar_periodic(inWin = vindauga, inShell = skall,
                                                         insDiv = insDiv,
                                                         inDims = c(xdim, ydim, dim(vindauga)[2], dim(skall)[2]))
    temp  <- matrix(data = temp, nrow = xdim, ncol = ydim)
    covariatesToCombine[[t]] <- temp[xlim[1]:xlim[2], ylim[1]:ylim[2]]
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdimShrunk, ydim = ydimShrunk)

  Z <- spatstat.geom::im(t(combinedCovariates), xrange = c(0, 1), yrange = c(0, numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, sDiv, r, h)

  return(Z)
}


#' Deprecated: function for creating beta squared covariates for the desired indices for DNS data
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param h lag
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @keywords internal
#' @export
combine.covariates.betasq.dns.shrunk <- function(sDiv, indices, h, xlim, ylim){
  warn.stacked.deprecated("combine.covariates.betasq.dns.shrunk")
  numberOfPictures <- length(indices)

  xdimShrunk <- xlim[2] - xlim[1] + 1
  ydimShrunk <- ylim[2] - ylim[1] + 1

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (t in 1:numberOfPictures) {
    insDiv <- sDiv[indices[t]+h, xlim[1]:xlim[2], ylim[1]:ylim[2]]
    covariatesToCombine[[t]] <- sum(insDiv^2) * matrix(data = 1, nrow = xdimShrunk, ncol = ydimShrunk)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdimShrunk, ydim = ydimShrunk)
  Z <- spatstat.geom::im(t(combinedCovariates), xrange = c(0, 1), yrange = c(0, numberOfPictures))

  rm(covariatesToCombine, xdimShrunk, ydimShrunk, numberOfPictures, sDiv, h)
  return(Z)
}

#' Deprecated: function for making a new AIC matrix given DNS data for a given radius and lag
#' @param r_and_h_untransformed r_and_h_combined in single number
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param sDiv horisontal divergence
#' @param indices which indices to use
#' @param radiiForPrint which radii to print
#' @param startLag First lag, used to decide when to print the radius
#' @param radiiModulo Parameter for how to decode the r_and_h_untransformed to r and h
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#' @param transformation Should the covariate be transformed, if value is "log" then a log transformation is applied, if value is "both", then the output will be both the transformation for the log and for identity
#'
#' @keywords internal
#' @export
aic.given.radius.and.lag.shrunk <- function(r_and_h_untransformed, vortices.pp, scars.pp, sDiv,
                                     indices, radiiForPrint, startLag, radiiModulo = 1000L,
                                     xlim, ylim, transformation = "identity"){
  warn.stacked.deprecated("aic.given.radius.and.lag.shrunk")
  r <- r_and_h_untransformed %% radiiModulo
  h <- (r_and_h_untransformed - r)/radiiModulo

  if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  Z <- combine.covariates.dns.shrunk(sDiv = sDiv, indices = indices, r = r, h = h, xlim = xlim, ylim = ylim)
  if(transformation == "log"){Z <- log(Z)}
  if(transformation == "both"){
    vort.fit.id <- spatstat.model::ppm(vortices.pp ~ Z)
    scar.fit.id <- spatstat.model::ppm(scars.pp ~ Z)

    vort.fit.log <- spatstat.model::ppm(vortices.pp ~ log(Z))
    scar.fit.log <- spatstat.model::ppm(scars.pp ~ log(Z))

    return(c(stats::AIC(vort.fit.id), stats::AIC(scar.fit.id), stats::AIC(vort.fit.log), stats::AIC(scar.fit.log)))
  }
  vort.fit <- spatstat.model::ppm(vortices.pp ~ Z)
  scar.fit <- spatstat.model::ppm(scars.pp ~ Z)

  rm(Z, r, h)
  return(c(stats::AIC(vort.fit), stats::AIC(scar.fit)))
}


#' Deprecated: function for making a new AIC matrix given DNS data
#' @param radiar Vector with all radii for local variance
#' @param lags Vector with all lag values to use
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param sDiv horisontal divergence
#' @param indices which indices to use
#' @param radiiModulo parameter for encoding lag and radius into single matrix
#' @param radiiForPrint Boolean variable to indicate whether to print progress or not
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#' @param transformation Should the covariate be transformed, if value is "log" then a log transformation is applied, if value is "both", then the output will be both the transformation for the log and for identity
#'
#' @keywords internal
#' @export
make.new.aic.matrix.dns.shrunk <- function(radiar, lags, vortices.pp, scars.pp, sDiv,
                                    indices, radiiModulo = 1000L, radiiForPrint = c(1),
                                    xlim, ylim, transformation = "identity"){
  warn.stacked.deprecated("make.new.aic.matrix.dns.shrunk")
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)
  aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = aic.given.radius.and.lag.shrunk, vortices.pp = vortices.pp,
                       scars.pp = scars.pp, sDiv = sDiv, indices = indices, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, xlim = xlim, ylim = ylim, transformation = transformation)

  aicMatrix.vortices <- aicMatrices[1, , ]
  aicMatrix.scars <- aicMatrices[2, , ]
  if(transformation == "both"){
    aicMatrix.vortices.log <- aicMatrices[3, , ]
    aicMatrix.scars.log <- aicMatrices[4, , ]
    rm(aicMatrices)
    return(list("aicMatrix.vortices.id" = aicMatrix.vortices, "aicMatrix.scars.id" = aicMatrix.scars))
  }

  rm(aicMatrices)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}
