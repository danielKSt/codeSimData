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

#' Function for making a new AIC matrix given experimental data
#' @param radiar Vector with all radii for local variance
#' @param hmax Maximal lag
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param hDiv horisontal divergence
#' @param indices which indices to use
#' @param ensembles vector of ensembles to use
#' @param printProgress Boolean variable to indicate whether to print progress or not
#' @param skim number of pixles to skim of the edge
#'
#' @export
make.new.aic.matrix.experimental <- function(radiar, hmax, vortices.pp, scars.pp, hDiv,
                                indices, ensembles, printProgress = FALSE, skim = 0){
  aicMatrix.scars <- matrix(0, nrow = 2*hmax+1, ncol = length(radiar))
  aicMatrix.vortices <- matrix(0, nrow = 2*hmax+1, ncol = length(radiar))

  for (r in radiar) {
    for (h in c(-hmax:hmax)) {
      if(printProgress){print(c(r, h))}
      Z <- combine.covariates.experimental(hDiv = hDiv, indices = indices, ensembles = ensembles, r = r, h = h, skim = skim)
      a <- spatstat.model::ppm(vortices.pp ~ Z)
      aicMatrix.vortices[h+hmax+1, r] <- stats::AIC(a)
      b <- spatstat.model::ppm(scars.pp ~ Z)
      aicMatrix.scars[h+hmax+1, r] <- stats::AIC(b)
    }
  }

  rm(a, b, Z, r, h)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}


#' Function for expanding an AIC matrix given experimental data
#' @param gamleRadiar Vector with all radii for local variance in the ole AIC matrices
#' @param nyeRadiar Vector with radii for local variance in expansion
#' @param hmax Maximal lag
#' @param gamleMatriser list with the AIC matrices that should be expanded
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param hDiv horisontal divergence
#' @param indices which indices to use
#' @param ensembles vector of ensembles to use
#' @param printProgress Boolean variable to indicate whether to print progress or not
#' @param skim number of pixles to skim of the edge
#'
#' @export
expand.aic.matrix.experimental <- function(gamleRadiar, nyeRadiar, hmax, gamleMatriser, vortices.pp, scars.pp, hDiv,
                              indices, ensembles, printProgress = TRUE, skim = 0){

  aicMatrix.scars <- matrix(0, nrow = 2*hmax+1, ncol = length(gamleRadiar)+length(nyeRadiar))
  aicMatrix.vortices <- matrix(0, nrow = 2*hmax+1, ncol = length(gamleRadiar)+length(nyeRadiar))

  aicMatrix.scars[1:(2*hmax+1), gamleRadiar] <- gamleMatriser$aicMatrix.scars
  aicMatrix.vortices[1:(2*hmax+1), gamleRadiar] <- gamleMatriser$aicMatrix.vortices

  for (r in nyeRadiar) {
    for (h in c(-hmax:hmax)) {
      if(printProgress){print(c(r, h))}
      Z <- combine.covariates.experimental(hDiv = hDiv, indices = indices, ensembles = ensembles, r = r, h = h, skim = skim)
      a <- spatstat.model::ppm(vortices.pp ~ Z)
      aicMatrix.vortices[h+hmax+1, r] <- stats::AIC(a)
      b <- spatstat.model::ppm(scars.pp ~ Z)
      aicMatrix.scars[h+hmax+1, r] <- stats::AIC(b)
    }
  }
  rm(a, b, Z, r, h)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
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
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param r radius to use in calculating local variance
#' @param h lag
#'
#' @export
combine.covariates.dns <- function(sDiv, indices, r, h){
  numberOfPictures <- length(indices)

  xdim <- dim(sDiv)[2]
  ydim <- dim(sDiv)[3]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  vindauga <- t(finn.vindauga.sirkel(radius = r))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for (t in 1:numberOfPictures) {
    insDiv <- sDiv[indices[t]+h, , ]
    covariatesToCombine[[t]] <- calculateLocVar_periodic(inWin = vindauga, inShell = skall,
                                                         insDiv = insDiv,
                                                         inDims = c(xdim, ydim, dim(vindauga)[2], dim(skall)[2]))
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim, ydim = ydim)

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdim),
                         seq(from = 0, to = numberOfPictures, length = ydim*numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, sDiv, r, h)

  return(Z)
}


#' Function for creating beta squared covariates for the desired indices for DNS data
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param h lag
#'
#' @export
combine.covariates.betasq.dns <- function(sDiv, indices, h){
  numberOfPictures <- length(indices)

  xdim <- dim(sDiv)[2]
  ydim <- dim(sDiv)[3]

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (t in 1:numberOfPictures) {
    insDiv <- sDiv[indices[t]+h, , ]
    covariatesToCombine[[t]] <- sum(insDiv^2) * matrix(data = 1, nrow = xdim, ncol = ydim)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdim, ydim = ydim)

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdim),
                         seq(from = 0, to = numberOfPictures, length = ydim*numberOfPictures))
  rm(covariatesToCombine, xdim, ydim, numberOfPictures, sDiv, h)

  return(Z)
}


#' Function for making a new AIC matrix given DNS data
#' @param radiar Vector with all radii for local variance
#' @param lags Vector with all lag values to use
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param sDiv horisontal divergence
#' @param indices which indices to use
#' @param printProgress Boolean variable to indicate whether to print progress or not
#'
#' @export
make.new.aic.matrix.dns <- function(radiar, lags, vortices.pp, scars.pp, sDiv,
                                             indices, printProgress = FALSE){
  aicMatrix.scars <- matrix(0, nrow = length(lags), ncol = length(radiar))
  aicMatrix.vortices <- matrix(0, nrow = length(lags), ncol = length(radiar))

  for (r in radiar) {
    for (h in lags) {
      if(printProgress){print(c(r, h))}
      Z <- combine.covariates.dns(sDiv = sDiv, indices = indices, r = r, h = h)
      a <- spatstat.model::ppm(vortices.pp ~ Z)
      aicMatrix.vortices[which(lags == h), r] <- stats::AIC(a)
      b <- spatstat.model::ppm(scars.pp ~ Z)
      aicMatrix.scars[which(lags == h), r] <- stats::AIC(b)
    }
  }

  rm(a, b, Z, r, h)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}

#' Function for expanding an AIC matrix given dns data
#' @param gamleRadiar Vector with all radii for local variance in the ole AIC matrices
#' @param nyeRadiar Vector with radii for local variance in expansion
#' @param lags Vector with all lag values to use
#' @param gamleMatriser list with the AIC matrices that should be expanded
#' @param vortices.pp votrices as spatstat point process
#' @param scars.pp scars as spatstat point process
#' @param sDiv horisontal divergence
#' @param indices which indices to use
#' @param printProgress Boolean variable to indicate whether to print progress or not
#'
#' @export
expand.aic.matrix.dns <- function(gamleRadiar, nyeRadiar, lags, gamleMatriser, vortices.pp, scars.pp, sDiv,
                                           indices, printProgress = TRUE){

  aicMatrix.scars <- matrix(0, nrow = length(lags), ncol = length(gamleRadiar)+length(nyeRadiar))
  aicMatrix.vortices <- matrix(0, nrow = length(lags), ncol = length(gamleRadiar)+length(nyeRadiar))

  aicMatrix.scars[1:length(lags), gamleRadiar] <- gamleMatriser$aicMatrix.scars
  aicMatrix.vortices[1:length(lags), gamleRadiar] <- gamleMatriser$aicMatrix.vortices

  for (r in nyeRadiar) {
    for (h in lags) {
      if(printProgress){print(c(r, h))}
      Z <- combine.covariates.dns(sDiv = sDiv, indices = indices, r = r, h = h)
      a <- spatstat.model::ppm(vortices.pp ~ Z)
      aicMatrix.vortices[which(lags == h), r] <- stats::AIC(a)
      b <- spatstat.model::ppm(scars.pp ~ Z)
      aicMatrix.scars[which(lags == h), r] <- stats::AIC(b)
    }
  }
  rm(a, b, Z, r, h)
  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}
