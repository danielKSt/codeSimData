#' Functions used for going from the complete mapped vortex data, to coordinate data for vortex centroids

#' @param pixelCoords coordinates of pixels in current dimension
#' @param modNum Total number of pixels in dimension
myMean.pixels.crossing.border <- function(pixelCoords, modNum){
  pixelCoords <- pixelCoords*(pixelCoords < floor(modNum/2)) +
    (pixelCoords - modNum)*(pixelCoords >= floor(modNum/2))
  return(c(mean(pixelCoords), max(pixelCoords)-min(pixelCoords)))
}

#' Function used to find vortex properties given all the pixels in the vortex
#' @param vortexPixels Data Frame with coordinates of all the pixels in the
#' @param boxDimensions vector of dimension of observation box
#'
#' @returns Vector of vortex properties
#' @export
explore.vortex <- function(vortexPixels, boxDimensions){
  if(any(vortexPixels$x==1) && any(vortexPixels$x==boxDimensions[1])){
    xInfo <- myMean.pixels.crossing.border(vortexPixels$x, boxDimensions[1])
  }
  else {xInfo <- c(mean(vortexPixels$x), max(vortexPixels$x)-min(vortexPixels$x))}

  if(any(vortexPixels$y==1) && any(vortexPixels$y==boxDimensions[2])){
    yInfo <- myMean.pixels.crossing.border(vortexPixels$y, boxDimensions[2])
  }
  else {yInfo <- c(mean(vortexPixels$y), max(vortexPixels$y)-min(vortexPixels$y))}
  return(c(xInfo[1]/(boxDimensions[1]+1),
           yInfo[1]/(boxDimensions[2]+1),
           dim(vortexPixels)[1],
           xInfo[2],
           yInfo[2]))
}


#' Function to go from complete mapped vortex pixel data to the point process format of vortex data
#' @param vortexData Matrix of zeros and ones, where the ones are the pixels with a vortex on them
#'
#' @returns Dataframe with vortex data
#' @export
find.vortices.specific.time <- function(vortexData){
  numberOfVortices <- 0
  for (x in 1:dim(vortexData)[1]) {
    for (y in 1:dim(vortexData)[2]) {
      if(vortexData[x,y]==1){
        x_start <- x
        y_start <- y
        numberOfVortices <- numberOfVortices + 1
        remaining <- list(c(x_start, y_start))
        remainingCount <- 1
        vortexData[x_start, y_start] <- 2
        pixels <- data.frame(x = x_start, y = y_start)
        pixelsCount <- 1
        while(remainingCount > 0){
          x <- remaining[[remainingCount]][1]
          y <- remaining[[remainingCount]][2]
          remaining[[remainingCount]] <- NULL
          remainingCount <- remainingCount - 1
          kj_list <- list(c(-1, -1), c(-1, 0), c(-1, 1), c(0, -1), c(0, 1), c(1, -1), c(1, 0), c(1, 1))
          for (kj in kj_list) {
            xNew <- (x + kj[1])%%dim(vortexData)[1]
            xNew <- xNew + dim(vortexData)[1]*(xNew==0)
            yNew <- (y + kj[2])%%dim(vortexData)[2]
            yNew <- yNew + dim(vortexData)[2]*(yNew==0)
            if(vortexData[xNew, yNew] == 1){
              remainingCount <- remainingCount + 1
              remaining[[remainingCount]] <- c(xNew, yNew)
              vortexData[xNew, yNew] <- 2
              pixels[nrow(pixels)+1,] <- c(xNew, yNew)
              pixelsCount <- pixelsCount + 1
            }
          }
        }
        vortexSummary <- explore.vortex(pixels, c(dim(vortexData)[1], dim(vortexData)[2]))
        if(numberOfVortices == 1){
          vortices <- data.frame(x = vortexSummary[1], y = vortexSummary[2], noPixels = vortexSummary[3],
                                 xStretch = vortexSummary[4], yStretch = vortexSummary[5])
        } else {
          vortices[nrow(vortices)+1,] <- vortexSummary
        }
      }
    }
  }
  if(numberOfVortices > 0){
    return(vortices)
  } else {
    return(0)
  }
}




#' Function to create a complete timeseries of all the vortex point process data
#'
#' @param vortexData T by M by N binary matrix, where a one at (t, m, n) indicates there is a vortex at pixel (m,n) at timestep t.
#'
#' @returns List of dataframes of vortex data for each timestep
#' @export
find.vortices.over.time <- function(vortexData){
  endTimestep <- dim(vortexData)[1]
  vortices <- vector(mode = "list", length = endTimestep)
  for(i in 1:endTimestep) {
    vortices[[i]] <- find.vortices.specific.time(vortexData = vortexData[i, , ])
  }

  return(vortices)
}
