#' Functions used for going from the complete mapped vortex data, to coordinate data for vortex centroids

#' @param pixelCoords coordinates of pixels in current dimension
#' @param modNum Total number of pixels in dimension
myMean.pixels.crossing.border <- function(pixelCoords, modNum){
  pixelCoords <- pixelCoords*(pixelCoords < floor(modNum/2)) +
    (pixelCoords - modNum)*(pixelCoords >= floor(modNum/2))
  return(c(mean(pixelCoords)%%modNum, max(pixelCoords)-min(pixelCoords)))
}

#' Function for mapping a (possibly fractional) pixel index to a coordinate in the unit square
#'
#' Pixel k (k = 1, ..., N) covers the interval ((k-1)/N, k/N), so its centre is at (k-0.5)/N.
#' Index 0 is the same pixel as index N, which is why the result is taken modulo 1.
#' This is the same convention spatstat uses for an image made with `xrange = c(0, 1)`.
#'
#' @param pixelIndex Pixel index, or mean pixel index, in the current dimension
#' @param modNum Total number of pixels in the dimension
#'
#' @returns Coordinate in the interval [0, 1)
pixel.to.unit.coord <- function(pixelIndex, modNum){
  return(((pixelIndex - 0.5)/modNum) %% 1.0)
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
  return(c(pixel.to.unit.coord(xInfo[1], boxDimensions[1]),
           pixel.to.unit.coord(yInfo[1], boxDimensions[2]),
           dim(vortexPixels)[1],
           xInfo[2],
           yInfo[2]))
}


#' Function to go from complete mapped vortex pixel data to the point process format of vortex data
#' @param vortexData Matrix of zeros and ones, where the ones are the pixels with a vortex on them
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns Dataframe with vortex data
#' @export
find.vortices.specific.time <- function(vortexData, periodic = TRUE){
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
            if(periodic){
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
            } else {
              inWindow <- (x + kj[1] <= dim(vortexData)[1]) && (y + kj[2] <= dim(vortexData)[2]) && (y + kj[2] > 0) && (x + kj[1] > 0)
              if(inWindow){
                xNew <- x + kj[1]
                yNew <- y + kj[2]
                if(vortexData[xNew, yNew] == 1){
                  remainingCount <- remainingCount + 1
                  remaining[[remainingCount]] <- c(xNew, yNew)
                  vortexData[xNew, yNew] <- 2
                  pixels[nrow(pixels)+1,] <- c(xNew, yNew)
                  pixelsCount <- pixelsCount + 1
                }
              }
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
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns List of dataframes of vortex data for each timestep
#' @export
find.vortices.over.time <- function(vortexData, periodic = TRUE){
  endTimestep <- dim(vortexData)[1]
  vortices <- vector(mode = "list", length = endTimestep)
  for(i in 1:endTimestep) {
    vortices[[i]] <- find.vortices.specific.time(vortexData = vortexData[i, , ], periodic = periodic)
  }

  attr(vortices, "pixelConvention") <- "centre"
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
  # Undo the shift so the centre is given as an index of the original pixel grid
  if(x_shift){
    centre_x <- centre_x - floor(dims[1]/2)
  }
  if(y_shift){
    centre_y <- centre_y - floor(dims[2]/2)
  }
  centre_x <- pixel.to.unit.coord(centre_x, dims[1])
  centre_y <- pixel.to.unit.coord(centre_y, dims[2])
  return(c(centre_x, centre_y, noPixels, xStretch, yStretch))
}


#' Function to go from complete mapped scar pixel data to the point process format of scar data
#' @param scarData Matrix of zeros and ones, where the ones are the pixels with a scar on them
#' @param np numpy from package reticulate
#' @param skimage skimage
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns Dataframe with scar data
#' @export
find.scars.specific.time <- function(scarData, np, skimage, periodic = TRUE){
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
            if(periodic){
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
            } else {
              inWindow <- (x + kj[1] <= dim(scarData)[1]) && (y + kj[2] <= dim(scarData)[2]) && (y + kj[2] > 0) && (x + kj[1] > 0)
              if(inWindow){
                xNew <- x + kj[1]
                yNew <- y + kj[2]
                if(scarData[xNew, yNew] == 1){
                  remainingCount <- remainingCount + 1
                  remaining[[remainingCount]] <- c(xNew, yNew)
                  scarData[xNew, yNew] <- 2
                  pixels[nrow(pixels)+1,] <- c(xNew, yNew)
                  pixelsCount <- pixelsCount + 1
                }
              }
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
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns List of dataframes of scar data for each timestep
#' @export
find.scars.over.time <- function(scarData, periodic = TRUE){
  reticulate::py_require("scikit-image")
  skimage <- reticulate::import("skimage.morphology")
  np <- reticulate::import("numpy")
  endTimestep <- dim(scarData)[1]
  scars <- vector(mode = "list", length = endTimestep)
  for(i in 1:endTimestep) {
    scars[[i]] <- find.scars.specific.time(scarData = scarData[i, , ], np = np, skimage = skimage, periodic = periodic)
  }
  attr(scars, "pixelConvention") <- "centre"
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
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns Dataframe with scar data
#' @export
find.scars.centrelines.specific.time <- function(scarData, np = NULL, skimage = NULL, periodic = TRUE){
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
            if(periodic){
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
            } else {
              inWindow <- (x + kj[1] <= dim(scarData)[1]) && (y + kj[2] <= dim(scarData)[2]) && (y + kj[2] > 0) && (x + kj[1] > 0)
              if(inWindow){
                xNew <- x + kj[1]
                yNew <- y + kj[2]
                if(scarData[xNew, yNew] == 1){
                  remainingCount <- remainingCount + 1
                  remaining[[remainingCount]] <- c(xNew, yNew)
                  scarData[xNew, yNew] <- 2
                  pixels[nrow(pixels)+1,] <- c(xNew, yNew)
                  pixelsCount <- pixelsCount + 1
                }
              }
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


#' Function to move points found with the old coordinate convention to the pixel-centre convention
#'
#' Before October 2026, vortex centres were placed at m/(N+1) and scar centres at m/N (modulo 1),
#' where m is the pixel index of the centre. Both are now placed at (m-0.5)/N, the centre of the pixel
#' in an image with `xrange = c(0, 1)`. This function converts point sets found with the old convention,
#' so they do not have to be found again. For scars the conversion is exact when the number of pixels is even.
#' When it is odd, scars that crossed the periodic border already sat at the pixel centre under the old code,
#' so the conversion moves them half a pixel too far.
#'
#' Before May 2026, explore.scar returned pixel indices for scars that did not cross the periodic border,
#' and coordinates in the unit square for scars that did. Point sets made then were divided by the number of
#' pixels afterwards, so the scars crossing the border were divided twice and piled up in the first pixel
#' row and column. Set `repairSeam = TRUE` to undo this before converting.
#'
#' @param points_list List of points for each timestep, as output from find.vortices.over.time or find.scars.over.time.
#' Can also be a list of such lists (one per ensemble), or the dataframe for a single timestep.
#' @param dims Vector with the number of pixels in the x and y direction of the field the points were found in
#' @param type "vortices" if the points were found with find.vortices.over.time, "scars" if found with find.scars.over.time
#' @param repairSeam Set to TRUE for scars where the coordinates of scars crossing the border were divided by the
#' number of pixels twice. Only used when type is "scars".
#'
#' @returns The points in the same format, with x and y moved to the pixel-centre convention
#' @export
convert.points.to.pixel.centres <- function(points_list, dims, type, repairSeam = FALSE){
  if(identical(attr(points_list, "pixelConvention"), "centre")){
    stop("These points already use the pixel-centre convention.")
  }
  if(!(type %in% c("vortices", "scars"))){
    stop("type must be either \"vortices\" or \"scars\".")
  }
  if(repairSeam && type != "scars"){
    stop("repairSeam can only be used when type is \"scars\".")
  }
  # Applies fun to the dataframe of every timestep, skipping empty timesteps
  apply.to.points <- function(points, fun){
    if(is.data.frame(points)){
      return(fun(points))
    }
    for(i in seq_along(points)){
      if(!is.numeric(points[[i]])){
        points[[i]] <- apply.to.points(points[[i]], fun)
      }
    }
    return(points)
  }

  if(repairSeam){
    # A scar divided twice lies in [0, 1/N), where no correctly scaled scar can be except at 0.
    # Values exactly 1 come from scars in pixel N, which is the same as 0 on the torus.
    repair.coord <- function(coord, modNum){
      return(ifelse(coord*modNum < 1, coord*modNum, coord %% 1.0))
    }
    points_list <- apply.to.points(points_list, function(points){
      points$x <- repair.coord(points$x, dims[1])
      points$y <- repair.coord(points$y, dims[2])
      return(points)
    })
  }

  if(type == "scars"){
    # Old scar coordinates are whole pixel indices divided by the number of pixels
    noOffGrid <- 0
    apply.to.points(points_list, function(points){
      offGrid <- (abs(points$x*dims[1] - round(points$x*dims[1])) > 1e-8) |
        (abs(points$y*dims[2] - round(points$y*dims[2])) > 1e-8)
      noOffGrid <<- noOffGrid + sum(offGrid)
      return(points)
    })
    if(noOffGrid > 0){
      warning(paste(noOffGrid, "scar(s) are not on the old pixel grid (multiples of 1/N), so the conversion",
                    "will not put them on pixel centres. If they were divided by the number of pixels twice,",
                    "use repairSeam = TRUE."))
    }
  }

  convert.coord <- function(coord, modNum){
    if(type == "vortices"){
      return(pixel.to.unit.coord(coord*(modNum+1), modNum))
    } else {
      return((coord - 0.5/modNum) %% 1.0)
    }
  }
  points_list <- apply.to.points(points_list, function(points){
    points$x <- convert.coord(points$x, dims[1])
    points$y <- convert.coord(points$y, dims[2])
    return(points)
  })
  attr(points_list, "pixelConvention") <- "centre"
  return(points_list)
}
