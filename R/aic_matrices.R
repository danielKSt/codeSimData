#### Common functions experimental and DNS ####

#' Function for combining covariates from multiple snapshots into one single image
#' @param covariates list of all covariates to combine. Every element should be a matrix of the same dimension
#' @param xdim xdim of covariates
#' @param ydim ydim of covariates
#'
#' @export
combine.covariates.matrices <- function(covariates, xdim, ydim){
  covariate <- matrix(data = NA, nrow = xdim, ncol = ydim*length(covariates))
  for (i in 1:length(covariates)) {
    covariate[ , ((i-1)*ydim+1):(i*ydim)] <- covariates[[i]]
  }
  rm(covariates, xdim, ydim)
  invisible(covariate)
}

#### Experimental data ####

#' Function for scaling points to window when skimming edge pixles
#' @param points which points to scale
#' @param c_x how to scale x-axis
#' @param c_y how to scale y-axis
#'
scale.points.skim <- function(points, c_x, c_y){
  points$x <- c_x*(points$x-0.5)+0.5
  points$y <- c_y*(points$y-0.5)+0.5
  points <- points[points$x >= 0, ]
  points <- points[points$x <= 1, ]
  points <- points[points$y >= 0, ]
  points <- points[points$y <= 1, ]
  return(points)
}

#' Function for preparing the points, given a list of all points
#' @param indices indices to use points from
#' @param numberOfEnsembles vector with all ensembles to use
#' @param points The complete list of points at all times
#' @param xdim xdim of covariates used in scaling
#' @param ydim ydim of covariates used in scaling
#' @param skim how many pixles to skim off the edge?
#'
#' @export
prepare.points.experimental <- function(indices, numberOfEnsembles, points, xdim, ydim, skim = 0){
  c_x <- xdim/(xdim-2*skim)
  c_y <- ydim/(ydim-2*skim)

  ensemble <- 1
  t <- 1

  firstPointNotFound <- TRUE
  while(firstPointNotFound){
    if (!is.numeric(points[[ensemble]][[indices[t]]])) {
      temp <- rbind(points[[ensemble]][[indices[t]]])
      temp <- scale.points.skim(points = temp, c_x = c_x, c_y = c_y)
      if(dim(temp)[1] > 0){
        firstPointNotFound <- FALSE
        l <- dim(temp)[1]
        pointsForAnalysis <- rbind(temp+c(rep(0, l), rep(t-1+(ensemble-1)*length(indices), l)))
        break
      }
    }
    t <- t + 1
    if(t <- length(indices)){
      t <- 1
      ensemble <- ensemble + 1
      if(ensemble > numberOfEnsembles){
        break
      }
    }
  }

  timeOfFirstPoint <- t
  ensembleOfFirstPoint <- ensemble
  for (t in (timeOfFirstPoint+1):length(indices)) {
    if(!is.numeric(points[[ensembleOfFirstPoint]][[indices[t]]])){
      temp <- scale.points.skim(points = points[[ensembleOfFirstPoint]][[indices[t]]], c_x = c_x, c_y = c_y)
      if(dim(temp)[1] > 0){
        l <- dim(temp)[1]
        pointsForAnalysis <- rbind(pointsForAnalysis, temp+c(rep(0, l), rep(t-1+(ensembleOfFirstPoint-1)*length(indices), l)))
      }
    }
  }
  for (ensemble in (ensembleOfFirstPoint+1):numberOfEnsembles) {
    for (t in 1:length(indices)) {
      if(!is.numeric(points[[ensemble]][[indices[t]]])){
        temp <- scale.points.skim(points = points[[ensemble]][[indices[t]]], c_x = c_x, c_y = c_y)
        if(dim(temp)[1] > 0){
          l <- dim(temp)[1]
          pointsForAnalysis <- rbind(pointsForAnalysis, temp+c(rep(0, l), rep(t-1+(ensemble-1)*length(indices), l)))
        }
      }
    }
  }

  points.pp <- spatstat.geom::ppp(x = pointsForAnalysis$x, y = pointsForAnalysis$y, c(0,1), c(0,numberOfEnsembles*length(indices)))
  rm(pointsForAnalysis)

  return(points.pp)
}

#' Function for creating covariates for the desired indices
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param ensembles vector with all ensembles to use
#' @param r radius to use in calculating local variance
#' @param h lag
#' @param skim how many pixles to skim off the edge?
#'
#' @export
combine.covariates.experimental <- function(hDiv, indices, ensembles, r, h, skim = 0){
  numberOfPictures <- length(indices)*length(ensembles)

  xdim <- dim(hDiv)[3]
  ydim <- dim(hDiv)[4]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  vindauga <- t(finn.vindauga.sirkel(radius = r))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for (i in 1:numberOfPictures) {
    t <- ((i-1) %% length(indices)+1)
    ens <- (i-t)/length(indices) + 1
    insDiv <- hDiv[ens, indices[t]+h, skim:(xdim-skim), skim:(ydim-skim)]
    covariatesToCombine[[i]] <- calculateLocVar_physical(inWin = vindauga, inShell = skall,
                                                         insDiv = insDiv,
                                                         inDims = c(xdim-2*skim, ydim-2*skim, dim(vindauga)[2], dim(skall)[2]))
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim-2*skim, ydim = ydim-2*skim)

  Z <- spatstat.geom::im(t(combinedCovariates),
          seq(from = 0, to = 1,length = xdim-2*skim),
          seq(from = 0, to = numberOfPictures, length = (ydim-2*skim)*numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, hDiv)

  return(Z)
}

#' Function for creating the beta squared covariates for experimental data
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param ensembles vector with all ensembles to use
#' @param h lag
#' @param skim how many pixles to skim off the edge?
#'
#' @export
combine.covariates.betasq.experimental <- function(hDiv, indices, ensembles, h, skim = 0){
  numberOfPictures <- length(indices)*length(ensembles)

  xdim <- dim(hDiv)[3]
  ydim <- dim(hDiv)[4]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (i in 1:numberOfPictures) {
    t <- ((i-1) %% length(indices)+1)
    ens <- (i-t)/length(indices) + 1
    insDiv <- hDiv[ens, indices[t]+h, skim:(xdim-skim), skim:(ydim-skim)]
    covariatesToCombine[[i]] <- sum(insDiv^2)*matrix(data = 1, nrow = xdim-2*skim,
                                                     ncol = ydim-2*skim)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine,
                                                    xdim = xdim-2*skim, ydim = ydim-2*skim)

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdim-2*skim),
                         seq(from = 0, to = numberOfPictures, length = (ydim-2*skim)*numberOfPictures))
  rm(covariatesToCombine, xdim, ydim, numberOfPictures, hDiv)

  return(Z)
}


#' Function for making a new AIC matrix given DNS data for a given radius and lag
#' @param r_and_h_untransformed r_and_h_combined in single number
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param hDiv horisontal divergence
#' @param indices which indices to use
#' @param radiiForPrint which radii to print
#' @param startLag First lag, used to decide when to print the radius
#' @param radiiModulo Parameter for how to decode the r_and_h_untransformed to r and h
#' @param transformation Should the covariate be transformed, if value is "log" then a log transformation is applied
#' @param ensembles vector with all ensembles to use
#'
#' @export
aic.given.radius.and.lag.experimental <- function(r_and_h_untransformed, vortices.pp, scars.pp, hDiv, ensembles,
                                     indices, radiiForPrint, startLag, radiiModulo = 1000L, transformation = "identity"){
  r <- r_and_h_untransformed %% radiiModulo
  h <- (r_and_h_untransformed - r)/radiiModulo

  if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  Z <- combine.covariates.experimental(hDiv = hDiv, indices = indices, ensembles = ensembles, r = r, h = h)
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

#### DNS data ####

#' Function for preparing the points, given a list of all points
#' @param indices vector of indices/times to use points from
#' @param points The complete list of points at all times
#'
#' @export
prepare.points.dns <- function(indices, points){
  t <- 1
  firstPointNotFound <- TRUE
  while(firstPointNotFound){
    if (!is.numeric(points[[indices[t]]])) {
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
    if(!is.numeric(points[[indices[t]]])){
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

#' Function for creating local variance covariates for the desired indices for DNS data
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param r radius to use in calculating local variance
#' @param h lag
#'
#' @export
combine.covariates.dns <- function(hDiv, indices, r, h){
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

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdim),
                         seq(from = 0, to = numberOfPictures, length = ydim*numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, hDiv, r, h)

  return(Z)
}

#' Function for creating beta squared covariates for the desired indices for DNS data
#' @param hDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param h lag
#'
#' @export
combine.covariates.betasq.dns <- function(hDiv, indices, h){
  numberOfPictures <- length(indices)

  xdim <- dim(hDiv)[2]
  ydim <- dim(hDiv)[3]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (t in 1:numberOfPictures) {
    insDiv <- hDiv[indices[t]+h, , ]
    covariatesToCombine[[t]] <- sum(insDiv^2) * matrix(data = 1, nrow = xdim, ncol = ydim)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim, ydim = ydim)

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdim),
                         seq(from = 0, to = numberOfPictures, length = ydim*numberOfPictures))
  rm(covariatesToCombine, xdim, ydim, numberOfPictures, hDiv, h)

  return(Z)
}

#' Function for making a new AIC matrix given DNS data for a given radius and lag
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
#' @export
aic.given.radius.and.lag <- function(r_and_h_untransformed, vortices.pp, scars.pp, hDiv,
                                    indices, radiiForPrint, startLag, radiiModulo = 1000L, transformation = "identity"){
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

#' Function for making a matrix containing information both on the lag and radius to be used for each index in the aicMatrices
#' @param radiar Vector with all radii for local variance
#' @param lags Vector with all lag values to use
#' @param radiiModulo parameter for how to encode lag and radius
#'
#' @export
make.input.matrix <- function(radiar, lags, radiiModulo = 1000L){
  res <- matrix(data = as.integer(radiar), ncol = length(radiar), nrow = length(lags), byrow = TRUE)
  for (h_i in 1:length(lags)) {
    res[h_i, ] <- res[h_i, ] + radiiModulo*lags[h_i]
  }
  return(res)
}


#' Function for making a new AIC matrix given data
#' @param radiar Vector with all radii for local variance
#' @param lags Vector with all lag values to use
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param hDiv horisontal divergence
#' @param indices which indices to use
#' @param ensembles Vector of ensembles to make use of, set to NULL if the data is not split into ensembles
#' @param sim Boolean set to true if the data is simulated with periodic BCs
#' @param skim Not yet available
#' @param radiiModulo parameter for encoding lag and radius into single matrix
#' @param radiiForPrint Boolean variable to indicate whether to print progress or not
#' @param transformation Should the covariate be transformed, if value is "log" then a log transformation is applied
#'
#' @export
make.new.aic.matrix <- function(radiar, lags, vortices.pp, scars.pp, hDiv, indices, ensembles = NULL, sim = TRUE,
                                skim = 0, radiiForPrint, radiiModulo = 1000L, transformation = "identity"){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)

  if(sim){
    aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = aic.given.radius.and.lag, vortices.pp = vortices.pp,
                         scars.pp = scars.pp, hDiv = hDiv, indices = indices, radiiForPrint = radiiForPrint,
                         startLag = lags[1], radiiModulo = radiiModulo, transformation = transformation)
  } else {
    aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = aic.given.radius.and.lag.experimental, vortices.pp = vortices.pp,
                         scars.pp = scars.pp, hDiv = hDiv, indices = indices, radiiForPrint = radiiForPrint, ensembles = ensembles,
                         startLag = lags[1], radiiModulo = radiiModulo, transformation = transformation)
  }

  aicMatrix.vortices <- aicMatrices[1, , ]
  aicMatrix.scars <- aicMatrices[2, , ]
  if(transformation == "both"){
    aicMatrix.vortices.log <- aicMatrices[3, , ]
    aicMatrix.scars.log <- aicMatrices[4, , ]
    rm(aicMatrices)
    return(list("aicMatrix.vortices.id" = aicMatrix.vortices, "aicMatrix.scars.id" = aicMatrix.scars,
                "aicMatrix.vortices.log" = aicMatrix.vortices.log, "aicMatrix.scars.log" = aicMatrix.scars.log))
  }

  rm(aicMatrices)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}
