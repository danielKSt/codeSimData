# Preparing points for analysis ----

# load(file = "/Volumes/work/danielks/RE2500_WEinf/simulatedPoints.RDa")
# load(file = "/Volumes/work/danielks/RE2500_WEinf/scales.RDa")
#
# noFrames <- length(simulatedVortices)
# startTime <- as.integer(ceiling(timescale*2/3)+1)
# lags <- seq(from = 1L-startTime, to = startTime-1L, by = 3L)
# tidspunkt <- seq(from = startTime, to = noFrames-max(lags), by = as.integer(ceiling(timescale)))
# rm(startTime)
#
# xlim <- c(0.25, 0.75)
# ylim <- c(0.25, 0.75)
#
# points <- simulatedVortices[[1]]
# point.in.shrunk(points = points, xlim = xlim, ylim = ylim)
#
# points <- simulatedVortices[[3030]]
# point.in.shrunk(points = points, xlim = xlim, ylim = ylim)
#
# points <- simulatedVortices[[1]]
# find.points.in.shrunk(points = points, xlim = xlim, ylim = ylim)
#
# indices <- tidspunkt
# points <- simulatedVortices
#
# plot(prepare.points.dns.shrunk(indices = c (1, 100), points = simulatedScars, xlim = xlim, ylim = ylim))
#
# plot(x = simulatedScars[[1]]$x, y = simulatedScars[[1]]$y)
# lines(x = c(0.25, 0.25, 0.75, 0.75, 0.25), y = c(0.25, 0.75, 0.75, 0.25, 0.25))
# plot(x = simulatedScars[[100]]$x, y = simulatedScars[[100]]$y)
# lines(x = c(0.25, 0.25, 0.75, 0.75, 0.25), y = c(0.25, 0.75, 0.75, 0.25, 0.25))

#' Function for scaling points from a reduced observation window into a "full" observation window
#' @param points dataframe of points
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
scale.points.shrunk <- function(points, xlim, ylim){
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

#' Function for preparing the points, given a list of all points
#' @param indices vector of indices/times to use points from
#' @param points The complete list of points at all times
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @export
prepare.points.dns.shrunk <- function(indices, points, xlim, ylim){
  t <- 1
  firstPointNotFound <- TRUE
  while(firstPointNotFound){
    if (!is.numeric(points[[indices[t]]])) {
      if (point.in.shrunk(points = points[[indices[t]]], xlim = xlim, ylim = ylim)){
        temp <- points[[indices[t]]]
        temp <- temp[find.points.in.shrunk(points = temp, xlim = xlim, ylim = ylim), ]
        temp <- scale.points.shrunk(points = temp, xlim = xlim, ylim = ylim)
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
    if(!is.numeric(points[[indices[t]]])){
      if(point.in.shrunk(points = points[[indices[t]]], xlim = xlim, ylim = ylim)){
        temp <- points[[indices[t]]]
        temp <- temp[find.points.in.shrunk(points = temp, xlim = xlim, ylim = ylim), ]
        temp <- scale.points.shrunk(points = temp, xlim = xlim, ylim = ylim)
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

# Preparing covariates for analysis ----

# library("rhdf5")
# library(codeSimData)
# dataFolder <- "/Volumes/work/danielks/RE2500_WEinf/"
# load(file = "/Volumes/work/danielks/RE2500_WEinf/simulatedPoints.RDa")
# load(file = "/Volumes/work/danielks/RE2500_WEinf/scales.RDa")
#
# noFrames <- length(simulatedVortices)
# startTime <- as.integer(ceiling(timescale*2/3)+1)
# lags <- seq(from = 1L-startTime, to = startTime-1L, by = 3L)
# tidspunkt <- seq(from = startTime, to = noFrames-max(lags), by = as.integer(ceiling(timescale)))
# rm(startTime)
#
# hDiv <- h5read(file = "/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/true_sDiv.mat", name = 'horDiv')
#
# sDiv <- hDiv
# indices <- tidspunkt[1:4]
# indices <- c(1)
# r <- 10
# h <- 0
# xlim <- c(65, 256-64)
# ylim <- c(65, 256-64)
#
# a <- combine.covariates.dns.shrunk(sDiv = sDiv, indices = indices, r =r, h = h, xlim = xlim, ylim = ylim)
# plot(a)
#
# vindauga <- t(finn.vindauga.sirkel(radius = r))
# skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
# insDiv <- sDiv[1, , ]
# temp <- calculateLocVar_periodic(inWin = vindauga, inShell = skall,
#                                  insDiv = insDiv,
#                                  inDims = c(256, 256, dim(vindauga)[2], dim(skall)[2]))
# b <- matrix(data = temp, nrow = 256)
# b_red <- b[xlim[1]:xlim[2], ylim[1]:ylim[2]]
# image(b)
# image(b_red)
# sum(abs(b_red-t(a$v)))



#' Function for creating local variance covariates for the desired indices for DNS data
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param r radius to use in calculating local variance
#' @param h lag
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @export
combine.covariates.dns.shrunk <- function(sDiv, indices, r, h, xlim, ylim){
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

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdimShrunk),
                         seq(from = 0, to = numberOfPictures, length = ydimShrunk*numberOfPictures))
  rm(covariatesToCombine, vindauga, skall, xdim, ydim, numberOfPictures, sDiv, r, h)

  return(Z)
}


#' Function for creating beta squared covariates for the desired indices for DNS data
#' @param sDiv horrisontal divergence complete array
#' @param indices indices to use points from
#' @param h lag
#' @param xlim how much to shrink the window in the x direction
#' @param ylim how much to shrink the window in the y direction
#'
#' @export
combine.covariates.betasq.dns.shrunk <- function(sDiv, indices, h, xlim, ylim){
  numberOfPictures <- length(indices)

  xdimShrunk <- xlim[2] - xlim[1] + 1
  ydimShrunk <- ylim[2] - ylim[1] + 1

  covariatesToCombine <- vector(mode = "list", length = numberOfPictures)

  for (t in 1:numberOfPictures) {
    insDiv <- sDiv[indices[t]+h, xlim[1]:xlim[2], ylim[1]:ylim[2]]
    covariatesToCombine[[t]] <- sum(insDiv^2) * matrix(data = 1, nrow = xdimShrunk, ncol = ydimShrunk)
  }

  combinedCovariates <- combine.covariates.matrices(covariates = covariatesToCombine, xdim = xdimShrunk, ydim = ydimShrunk)

  Z <- spatstat.geom::im(t(combinedCovariates),
                         seq(from = 0, to = 1,length = xdimShrunk),
                         seq(from = 0, to = numberOfPictures, length = ydimShrunk*numberOfPictures))
  rm(covariatesToCombine, xdimShrunk, ydimShrunk, numberOfPictures, sDiv, h)

  return(Z)
}

#' Function for making a new AIC matrix given DNS data for a given radius and lag
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
#'
#' @export
aic.given.radius.and.lag.shrunk <- function(r_and_h_untransformed, vortices.pp, scars.pp, sDiv,
                                     indices, radiiForPrint, startLag, radiiModulo = 1000L,
                                     xlim, ylim){
  r <- r_and_h_untransformed %% radiiModulo
  h <- (r_and_h_untransformed - r)/radiiModulo

  if((r %in% radiiForPrint) && (h == startLag)) {print(r)}
  Z <- combine.covariates.dns.shrunk(sDiv = sDiv, indices = indices, r = r, h = h, xlim = xlim, ylim = ylim)
  vort.fit <- spatstat.model::ppm(vortices.pp ~ Z)
  scar.fit <- spatstat.model::ppm(scars.pp ~ Z)

  rm(Z, r, h)
  return(c(stats::AIC(vort.fit), stats::AIC(scar.fit)))
}


#' Function for making a new AIC matrix given DNS data
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
#'
#' @export
make.new.aic.matrix.dns.shrunk <- function(radiar, lags, vortices.pp, scars.pp, sDiv,
                                    indices, radiiModulo = 1000L, radiiForPrint = c(1),
                                    xlim, ylim){
  inputMatrix <- make.input.matrix(radiar = radiar, lags = lags, radiiModulo = radiiModulo)
  aicMatrices <- apply(X = inputMatrix, MARGIN = c(1, 2), FUN = aic.given.radius.and.lag.shrunk, vortices.pp = vortices.pp,
                       scars.pp = scars.pp, sDiv = sDiv, indices = indices, radiiForPrint = radiiForPrint,
                       startLag = lags[1], radiiModulo = radiiModulo, xlim = xlim, ylim = ylim)

  aicMatrix.vortices <- aicMatrices[1, , ]
  aicMatrix.scars <- aicMatrices[2, , ]
  rm(aicMatrices)

  return(list("aicMatrix.vortices" = aicMatrix.vortices, "aicMatrix.scars" = aicMatrix.scars))
}
