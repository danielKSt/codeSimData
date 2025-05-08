#' Function for distance between two points between 0 and 1 on the torus
#' @param x1 First number
#' @param x2 Second number
#'
dist.on.torus <- function(x1, x2){
  d <- abs(x1-x2)
  if(d <=0.5){
    return(d^2)
  } else {
    return((1-d)^2)
  }
}

#' Function for creating beta squared covariates for the desired indices for DNS data
#'
#' @param points.t0 dataframe with points to initialize
#'
initialize.point.tracking <- function(points.t0){
  pointsComplete <- vector(mode = "list", length = dim(points.t0)[1])
  living <- matrix(data = NA, nrow = length(pointsComplete), ncol = 2)
  for (i in 1:dim(points.t0)[1]) {
    pointsComplete[[i]] <- data.frame(t = 1, x = points.t0$x[i], y = points.t0$y[i], size = points.t0$noPixels[i], id = i)
    living[i, 1] <- i
    living[i, 2] <- i
  }
  return(list(pointsComplete = pointsComplete, living = living))
}

#' Function for finding the distance between all points at one timestep, and all points at the next timestep
#'
#' @param points.t1 dataframe with points at the first of the two timesteps
#' @param points.t2 dataframe with points at the second of the two timesteps
#' @param weight_dist parameter to indicate how to weight difference in distance between to points vs the difference in size of two points
#' @param moveOnTorus set to TRUE if points move on a torus
#'
closeness.of.points <- function(points.t1, points.t2, weight_dist = 256, moveOnTorus = TRUE){
  no_v1 <- dim(points.t1)[1]
  no_v2 <- dim(points.t2)[1]
  dist <- matrix(data = 0, nrow = no_v1, ncol = dim(points.t2)[1])
  for (i in 1:no_v1) {
    for (j in 1:no_v2) {
      if(moveOnTorus){
        dist[i,j] <- sqrt(dist.on.torus(points.t1$x[i],points.t2$x[j]) + dist.on.torus(points.t1$y[i],points.t2$y[j]))
      } else {
        dist[i,j] <- sqrt((points.t1$x[i] - points.t2$x[j])^2 + (points.t1$y[i] - points.t2$y[j])^2)
      }
      dist[i,j] <- weight_dist*dist[i,j] + abs(points.t1$noPixels[i]-points.t2$noPixels[j])
    }
  }
  return(dist)
}

#' Function for tracking points acros a single timestep after some of the initial precautions have been taken care of
#'
#' @param points.t1 dataframe with points at the first of the two timesteps
#' @param points.t2 dataframe with points at the second of the two timesteps
#' @param pointsComplete dataframe to put results into
#' @param living matrix with indices of the "living" vortices, first the indices in points.t1, then the related index in pointsComplete
#' @param t current timestep
#' @param closenessThreshold Parameter to control how far away vortices are allowed to be before they are judged to not be the same
#' @param weight_dist parameter to indicate how to weight difference in distance between to points vs the difference in size of two points
#' @param moveOnTorus set to TRUE if points move on a torus
#'
one.step <- function(points.t1, points.t2, pointsComplete, living, t, closenessThreshold = 20, weight_dist = 256, moveOnTorus = TRUE){
  no_v1 <- dim(points.t1)[1]
  no_v2 <- dim(points.t2)[1]

  dist <- closeness.of.points(points.t1 = points.t1, points.t2 = points.t2, moveOnTorus = moveOnTorus)
  minimal <- which(dist == min(dist) , arr.ind = TRUE)
  if(any(duplicated(minimal[ , 1]))){
    a  <- c(1:nrow(minimal))
    minimal <- matrix(data = minimal[-a[duplicated(minimal[, 1])], ], ncol = 2)
  }
  if(any(duplicated(minimal[ , 2]))){
    a  <- c(1:nrow(minimal))
    minimal <- matrix(data = minimal[-a[duplicated(minimal[, 2])], ], ncol = 2)
  }

  couples <- matrix(data = NA, nrow = min(no_v1, no_v2), ncol = 2)
  couples[1:dim(minimal)[1], ] <- minimal

  added.t1 <- minimal[, 1]
  added.t2 <- minimal[, 2]
  noAdded <- length(added.t1)

  while(noAdded < min(no_v1, no_v2) && min(dist[-added.t1, -added.t2]) < closenessThreshold){
    minimal <- which(dist == min(dist[-added.t1, -added.t2]), arr.ind = TRUE)
    newFirstConnector <- c(couples[1:noAdded , 1], minimal[ , 1])
    if(any(duplicated(newFirstConnector))){
      minimal <- matrix(data = minimal[-(which(duplicated(newFirstConnector))-noAdded), ], ncol = 2)
    }
    newSecondConnector <- c(couples[1:noAdded , 2], minimal[ , 2])
    if(any(duplicated(newSecondConnector))){
      minimal <- matrix(data = minimal[-(which(duplicated(newSecondConnector))-noAdded), ], ncol = 2)
    }
    if(any(duplicated(minimal[ , 1]))){
      a  <- c(1:nrow(minimal))
      minimal <- matrix(data = minimal[-a[duplicated(minimal[, 1])], ], ncol = 2)
    }
    if(any(duplicated(minimal[ , 2]))){
      a  <- c(1:nrow(minimal))
      minimal <- matrix(data = minimal[-a[duplicated(minimal[, 2])], ], ncol = 2)
    }

    couples[(noAdded+1):(noAdded + dim(minimal)[1]), ] <- minimal
    added.t1[(noAdded+1):(noAdded + dim(minimal)[1])] <- minimal[, 1]
    added.t2[(noAdded+1):(noAdded + dim(minimal)[1])] <- minimal[, 2]
    noAdded <- noAdded + dim(minimal)[1]
  }

  died <- setdiff(x = c(1:no_v1), y = couples[, 1])
  born <- setdiff(x = c(1:no_v2), y = couples[, 2])

  living_new <- matrix(data = 0, nrow = no_v2, ncol = 2)
  for (couple in 1:noAdded) {
    localID <- which(living[ , 1] == couples[couple, 1])
    pointsComplete[[living[localID, 2]]][nrow(pointsComplete[[living[localID, 2]]]) + 1,] <- c(t+1, points.t2[couples[couple, 2], 1:3], couples[couple, 2])
    living_new[couple, 1] <- couples[couple, 2]
    living_new[couple, 2] <- living[localID, 2]
  }

  id_start <- length(pointsComplete)
  if(length(born)>0){
    for(j in 1:length(born)){
      pointsComplete[[id_start+j]] <- data.frame(t = t+1, x = points.t2$x[born[j]], y = points.t2$y[born[j]],
                                                 size = points.t2$noPixels[born[j]], id = born[j])
      living_new[noAdded + j, 1] <- born[j]
      living_new[noAdded + j, 2] <- id_start + j
    }
    rm(j)
  }


  living <- living_new

  rm(born, died, living_new, id_start, points.t1, points.t2, couple, minimal, couples, no_v1, no_v2, localID)
  rm(added.t1, added.t2, noAdded, dist)

  return(list(pointsComplete = pointsComplete, living = living))
}

#' Function for tracking points across a single time step
#'
#' @param points.t1 dataframe with points at the first of the two timesteps
#' @param points.t2 dataframe with points at the second of the two timesteps
#' @param pointsComplete dataframe to put results into
#' @param living matrix with indices of the "living" vortices, first the indices in points.t1, then the related index in pointsComplete
#' @param t current timestep
#' @param closenessThreshold Parameter to control how far away vortices are allowed to be before they are judged to not be the same
#' @param weight_dist parameter to indicate how to weight difference in distance between to points vs the difference in size of two points
#' @param moveOnTorus set to TRUE if points move on a torus
#'
track.points.single.step <- function(points.t1, points.t2, pointsComplete, living, t, closenessThreshold = 20, weight_dist = 256, moveOnTorus = TRUE){
  no_v1 <- dim(points.t1)[1]
  no_v2 <- dim(points.t2)[1]
  if(is.numeric(points.t2)){
    living <- NULL
  } else if(is.numeric(points.t1)){
    id_start <- length(pointsComplete)
    living <- matrix(data = NA, nrow = no_v2, ncol = 2)
    for (j in 1:no_v2) {
      pointsComplete[[id_start+j]] <- data.frame(t = t+1, x = points.t2$x[j], y = points.t2$y[j],
                                                 size = points.t2$noPixels[j], id = id_start+j)
      living[j, 2] <- id_start + j
      living[j, 1] <- j
    }
  } else {
    res <- one.step(points.t1 = points.t1, points.t2 = points.t2, pointsComplete = pointsComplete, living = living,
                    t = t, closenessThreshold = closenessThreshold, weight_dist = weight_dist, moveOnTorus = moveOnTorus)
    pointsComplete <- res$pointsComplete
    living <- res$living
  }
  return(list(pointsComplete = pointsComplete, living = living))
}

#' Function for tracking points moving in some observation window
#' @param pointsUntracked A list of snapshots of points, each point on a snapshot must contain at least, x, y and noPixels
#' @param tEnd time to end
#' @param tStart time to start, if no input is given, this will be set to the length of pointsUntracked
#' @param closenessThreshold Parameter to control how far away vortices are allowed to be before they are judged to not be the same
#' @param weight_dist parameter to indicate how to weight difference in distance between to points vs the difference in size of two points
#' @param moveOnTorus set to TRUE if points move on a torus
#'
#' @export
track.points.over.time <- function(pointsUntracked, tEnd = NULL, tStart = 1, closenessThreshold = 20, weight_dist = 256, moveOnTorus = TRUE){
  t <- tStart
  while(is.numeric(pointsUntracked[[t]])){
    t <- t + 1
    if(t == tEnd){
      print("No vortices in the desired timeframe")
      return(NULL)
    }
  }
  tStart <- t
  res <- initialize.point.tracking(points.t0 = pointsUntracked[[tStart]])
  if(is.null(tEnd)){tEnd <- length(pointsUntracked)}
  for (t in tStart:(tEnd-1)) {
    res <- track.points.single.step(points.t1 = pointsUntracked[[t]], points.t2 = pointsUntracked[[t+1]],
                                    pointsComplete = res$pointsComplete, living = res$living, t = t,
                                    closenessThreshold = closenessThreshold, weight_dist = weight_dist, moveOnTorus = moveOnTorus)
  }
  return(res$pointsComplete)
}

