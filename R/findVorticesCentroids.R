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

# Scars below here:

#' Function used to find scar properties given all the pixels in the scar
#' @param pixels Data Frame with coordinates of all the pixels in the scar
#' @param dims vector of dimension of observation box
#' @param np numpy from package reticulate
#' @param skimage skimage
#'
#' @returns Vector of vortex properties
#' @export
explore.scar <- function(pixels, dims, np, skimage){
  if(is.null(np)){
    np <- reticulate::import("numpy")
  }
  if(is.null(skimage)){
    reticulate::py_require("scikit-image")
    skimage <- reticulate::import("skimage.morphology")
  }
  noPixels <- dim(pixels)[1]
  if((min(pixels$x) == 1) && (max(pixels$x) == dims[1])){
    pixels$x <- ((pixels$x + floor(dims[1]/2)) %% dims[1])
    x_shift <- TRUE
  } else {
    x_shift <- FALSE
  }
  if((min(pixels$y) == 1) && (max(pixels$y) == dims[2])){
    pixels$y <- ((pixels$y + floor(dims[2]/2)) %% dims[2])
    y_shift <- TRUE
  } else {
    y_shift <- FALSE
  }

  xStretch <- max(pixels$x) - min(pixels$x)
  yStretch <- max(pixels$y) - min(pixels$y)

  geom_centre_x <- mean(pixels$x)
  geom_centre_y <- mean(pixels$y)

  centreline <- skeletonize_object_py(df = pixels, np = np, skimage = skimage)
  minDistance <- dims[1]*2
  for (i in 1:nrow(centreline)) {
    dist <- sqrt((geom_centre_y-centreline$y[i])^2+(geom_centre_x-centreline$x[i])^2)
    if(dist < minDistance){
      minDistance <- dist
      centre_x <- centreline$x[i]
      centre_y <- centreline$y[i]
    }
  }
  if(x_shift){
    centre_x <- ((centre_x/dims[1])-0.5) %% 1.0
  }
  if(y_shift){
    centre_y <- ((centre_y/dims[2])-0.5) %% 1.0
  }
  return(c(centre_x, centre_y, noPixels, xStretch, yStretch))
}


#' Function to go from complete mapped scar pixel data to the point process format of scar data
#' @param scarData Matrix of zeros and ones, where the ones are the pixels with a scar on them
#' @param np numpy from package reticulate
#' @param skimage skimage
#'
#' @returns Dataframe with scar data
#' @export
find.scars.specific.time <- function(scarData, np, skimage){
  if(is.null(np)){
    np <- reticulate::import("numpy")
  }
  if(is.null(skimage)){
    reticulate::py_require("scikit-image")
    skimage <- reticulate::import("skimage.morphology")
  }
  numberOfScars <- 0
  for (x in 1:dim(scarData)[1]) {
    for (y in 1:dim(scarData)[2]) {
      if(scarData[x,y]==1){
        x_start <- x
        y_start <- y
        numberOfScars <- numberOfScars + 1
        remaining <- list(c(x_start, y_start))
        remainingCount <- 1
        scarData[x_start, y_start] <- 2
        pixels <- data.frame(x = x_start, y = y_start)
        pixelsCount <- 1
        while(remainingCount > 0){
          x <- remaining[[remainingCount]][1]
          y <- remaining[[remainingCount]][2]
          remaining[[remainingCount]] <- NULL
          remainingCount <- remainingCount - 1
          kj_list <- list(c(-1, -1), c(-1, 0), c(-1, 1), c(0, -1), c(0, 1), c(1, -1), c(1, 0), c(1, 1))
          for (kj in kj_list) {
            xNew <- (x + kj[1])%%dim(scarData)[1]
            xNew <- xNew + dim(scarData)[1]*(xNew==0)
            yNew <- (y + kj[2])%%dim(scarData)[2]
            yNew <- yNew + dim(scarData)[2]*(yNew==0)
            if(scarData[xNew, yNew] == 1){
              remainingCount <- remainingCount + 1
              remaining[[remainingCount]] <- c(xNew, yNew)
              scarData[xNew, yNew] <- 2
              pixels[nrow(pixels)+1,] <- c(xNew, yNew)
              pixelsCount <- pixelsCount + 1
            }
          }
        }
        scar <- explore.scar(pixels, c(dim(scarData)[1], dim(scarData)[2]), np = np, skimage = skimage)
        if(numberOfScars == 1){
          scars <- data.frame(x = scar[1], y = scar[2], noPixels = scar[3],
                              xStretch = scar[4], yStretch = scar[5])
        } else {
          scars[numberOfScars, ] <- scar
        }
      }
    }
  }
  if(numberOfScars > 0){
    return(scars)
  } else {
    return(0)
  }
}




#' Function to create a complete timeseries of all the scar point process data
#'
#' @param scarData T by M by N binary array, where a one at (t, m, n) indicates there is a scar at pixel (m,n) at timestep t.
#'
#' @returns List of dataframes of scar data for each timestep
#' @export
find.scars.over.time <- function(scarData){
  reticulate::py_require("scikit-image")
  skimage <- reticulate::import("skimage.morphology")
  np <- reticulate::import("numpy")
  endTimestep <- dim(scarData)[1]
  scars <- vector(mode = "list", length = endTimestep)
  for(i in 1:endTimestep) {
    scars[[i]] <- find.scars.specific.time(scarData = scarData[i, , ], np = np, skimage = skimage)
  }

  return(scars)
}


#' Skeletonize a scar given pixels
#'
#' @param df dataframe with the pixel coordinates of the scar
#' @param np numpy from package reticulate
#' @param skimage skimage
#'
#' @returns Centerline of the scar
#'
skeletonize_object_py <- function(df, np, skimage) {
  if(is.null(np)){
    np <- reticulate::import("numpy")
  }
  if(is.null(skimage)){
    reticulate::py_require("scikit-image")
    skimage <- reticulate::import("skimage.morphology")
  }
  # Create binary image matrix
  width <- max(df$x) + 1
  height <- max(df$y) + 1
  mat <- matrix(0L, nrow = height, ncol = width)
  for (i in 1:nrow(df)) {
    mat[df$y[i] + 1, df$x[i] + 1] <- 1L
  }

  # Convert to NumPy array
  np_mat <- np$array(t(mat))  # transpose to match Python's row-major order

  # Call skeletonize from scikit-image
  skeleton_np <- skimage$skeletonize(np_mat)

  # Convert back to R matrix
  skeleton_mat <- t(reticulate::py_to_r(skeleton_np))  # transpose back

  # Extract (x, y) coordinates where skeleton is True
  coords <- which(skeleton_mat == 1, arr.ind = TRUE)
  skeleton_df <- data.frame(
    x = coords[, "col"] - 1,
    y = coords[, "row"] - 1
  )

  return(skeleton_df)
}


#' Function to go from complete mapped scar pixel data to the point process format of scar data
#' @param scarData Matrix of zeros and ones, where the ones are the pixels with a scar on them
#' @param np numpy from package reticulate
#' @param skimage skimage
#'
#' @returns Dataframe with scar data
#' @export
find.scars.centrelines.specific.time <- function(scarData, np = NULL, skimage = NULL){
  if(is.null(np)){
    np <- reticulate::import("numpy")
  }
  if(is.null(skimage)){
    reticulate::py_require("scikit-image")
    skimage <- reticulate::import("skimage.morphology")
  }
  numberOfScars <- 0
  for (x in 1:dim(scarData)[1]) {
    for (y in 1:dim(scarData)[2]) {
      if(scarData[x,y]==1){
        x_start <- x
        y_start <- y
        numberOfScars <- numberOfScars + 1
        remaining <- list(c(x_start, y_start))
        remainingCount <- 1
        scarData[x_start, y_start] <- 2
        pixels <- data.frame(x = x_start, y = y_start)
        pixelsCount <- 1
        while(remainingCount > 0){
          x <- remaining[[remainingCount]][1]
          y <- remaining[[remainingCount]][2]
          remaining[[remainingCount]] <- NULL
          remainingCount <- remainingCount - 1
          kj_list <- list(c(-1, -1), c(-1, 0), c(-1, 1), c(0, -1), c(0, 1), c(1, -1), c(1, 0), c(1, 1))
          for (kj in kj_list) {
            xNew <- (x + kj[1])%%dim(scarData)[1]
            xNew <- xNew + dim(scarData)[1]*(xNew==0)
            yNew <- (y + kj[2])%%dim(scarData)[2]
            yNew <- yNew + dim(scarData)[2]*(yNew==0)
            if(scarData[xNew, yNew] == 1){
              remainingCount <- remainingCount + 1
              remaining[[remainingCount]] <- c(xNew, yNew)
              scarData[xNew, yNew] <- 2
              pixels[nrow(pixels)+1,] <- c(xNew, yNew)
              pixelsCount <- pixelsCount + 1
            }
          }
        }
        if(numberOfScars == 1){
          scars <- list(skeletonize_object_py(df = pixels, np = np, skimage = skimage))
        } else {
          scars[[numberOfScars]] <- skeletonize_object_py(df = pixels, np = np, skimage = skimage)
        }
      }
    }
  }
  if(numberOfScars > 0){
    return(scars)
  } else {
    return(0)
  }
}
