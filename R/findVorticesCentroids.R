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


#' Function for finding the connected structures in a binary mask
#' @description
#' Two pixels belong to the same structure if they are neighbours, including diagonal neighbours. With periodic
#' boundary conditions, pixels on opposite edges of the mask are neighbours too. The structures are ordered by their
#' first pixel when the mask is scanned with the x index in the outer loop.
#'
#' @param mask Matrix of zeros and ones, where the ones are the pixels with a structure on them
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns List with one data frame per structure, with the x and y indices of its pixels
find.connected.components <- function(mask, periodic = TRUE){
  nx <- dim(mask)[1]
  ny <- dim(mask)[2]
  unvisited <- mask == 1
  offsets <- rbind(c(-1, -1), c(-1, 0), c(-1, 1), c(0, -1), c(0, 1), c(1, -1), c(1, 0), c(1, 1))
  starts <- which(unvisited, arr.ind = TRUE)
  starts <- starts[order(starts[, 1], starts[, 2]), , drop = FALSE]
  components <- list()
  for(s in seq_len(nrow(starts))){
    if(!unvisited[starts[s, 1], starts[s, 2]]){
      next
    }
    unvisited[starts[s, 1], starts[s, 2]] <- FALSE
    # Breadth-first search, with the pixels stored as linear indices into the mask
    queue <- starts[s, 1] + nx*(starts[s, 2] - 1)
    head <- 1
    while(head <= length(queue)){
      x <- (queue[head] - 1) %% nx + 1
      y <- (queue[head] - 1) %/% nx + 1
      head <- head + 1
      xNew <- x + offsets[, 1]
      yNew <- y + offsets[, 2]
      if(periodic){
        xNew <- (xNew - 1) %% nx + 1
        yNew <- (yNew - 1) %% ny + 1
      } else {
        inWindow <- (xNew >= 1) & (xNew <= nx) & (yNew >= 1) & (yNew <= ny)
        xNew <- xNew[inWindow]
        yNew <- yNew[inWindow]
      }
      # unique() matters for masks one or two pixels wide, where several neighbours wrap to the same pixel
      neighbours <- unique(xNew + nx*(yNew - 1))
      neighbours <- neighbours[unvisited[neighbours]]
      unvisited[neighbours] <- FALSE
      queue <- c(queue, neighbours)
    }
    components[[length(components) + 1]] <- data.frame(x = (queue - 1) %% nx + 1, y = (queue - 1) %/% nx + 1)
  }
  return(components)
}

#' Function for collecting the properties of each structure in a data frame
#' @param summaries List with one vector per structure, as returned by explore.vortex or explore.scar
#'
#' @returns Data frame with one row per structure, and zero rows if there are no structures
structures.to.data.frame <- function(summaries){
  res <- matrix(data = as.numeric(unlist(summaries)), ncol = 5, byrow = TRUE,
                dimnames = list(NULL, c("x", "y", "noPixels", "xStretch", "yStretch")))
  return(as.data.frame(res))
}

#' Function for checking if a timestep has no structures
#' @description
#' Timesteps without structures are data frames with zero rows. Point sets found with older versions of the package
#' store them as the number 0 instead, so both count as empty.
#'
#' @param points Data frame with the structures of one timestep, or 0
#'
#' @returns TRUE if there are no structures
is.empty.snapshot <- function(points){
  return(is.numeric(points) || (is.data.frame(points) && nrow(points) == 0))
}

#' Function to go from complete mapped vortex pixel data to the point process format of vortex data
#' @param vortexData Matrix of zeros and ones, where the ones are the pixels with a vortex on them
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns Dataframe with vortex data, with zero rows if there are no vortices
#' @export
find.vortices.specific.time <- function(vortexData, periodic = TRUE){
  components <- find.connected.components(mask = vortexData, periodic = periodic)
  return(structures.to.data.frame(lapply(components, explore.vortex, boxDimensions = dim(vortexData))))
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

#' Function for getting numpy and skimage.morphology from Python
#' @description
#' Imports the modules through reticulate, unless they are given. scikit-image is required before the first import, so
#' that it is part of the Python environment reticulate sets up.
#'
#' @param np numpy from package reticulate, or NULL to import it
#' @param skimage skimage.morphology from package reticulate, or NULL to import it
#'
#' @returns List with np and skimage
python.modules <- function(np = NULL, skimage = NULL){
  if(is.null(skimage)){
    reticulate::py_require("scikit-image")
    skimage <- reticulate::import("skimage.morphology")
  }
  if(is.null(np)){
    np <- reticulate::import("numpy")
  }
  return(list(np = np, skimage = skimage))
}

#' Function used to find scar properties given all the pixels in the scar
#' @param pixels Data Frame with coordinates of all the pixels in the scar
#' @param dims vector of dimension of observation box
#' @param np numpy from package reticulate, or NULL to import it
#' @param skimage skimage.morphology from package reticulate, or NULL to import it
#'
#' @returns Vector of vortex properties
#' @export
explore.scar <- function(pixels, dims, np = NULL, skimage = NULL){
  modules <- python.modules(np = np, skimage = skimage)
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

  centreline <- skeletonize_object_py(df = pixels, np = modules$np, skimage = modules$skimage)
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
#' @param np numpy from package reticulate, or NULL to import it
#' @param skimage skimage.morphology from package reticulate, or NULL to import it
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns Dataframe with scar data, with zero rows if there are no scars
#' @export
find.scars.specific.time <- function(scarData, np = NULL, skimage = NULL, periodic = TRUE){
  modules <- python.modules(np = np, skimage = skimage)
  components <- find.connected.components(mask = scarData, periodic = periodic)
  return(structures.to.data.frame(lapply(components, explore.scar, dims = dim(scarData), np = modules$np,
                                         skimage = modules$skimage)))
}




#' Function to create a complete timeseries of all the scar point process data
#'
#' @param scarData T by M by N binary array, where a one at (t, m, n) indicates there is a scar at pixel (m,n) at timestep t.
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns List of dataframes of scar data for each timestep
#' @export
find.scars.over.time <- function(scarData, periodic = TRUE){
  modules <- python.modules()
  endTimestep <- dim(scarData)[1]
  scars <- vector(mode = "list", length = endTimestep)
  for(i in 1:endTimestep) {
    scars[[i]] <- find.scars.specific.time(scarData = scarData[i, , ], np = modules$np, skimage = modules$skimage,
                                           periodic = periodic)
  }
  attr(scars, "pixelConvention") <- "centre"
  return(scars)
}


#' Skeletonize a scar given pixels
#'
#' @param df dataframe with the pixel coordinates of the scar
#' @param np numpy from package reticulate, or NULL to import it
#' @param skimage skimage.morphology from package reticulate, or NULL to import it
#'
#' @returns Centerline of the scar
#'
skeletonize_object_py <- function(df, np = NULL, skimage = NULL) {
  modules <- python.modules(np = np, skimage = skimage)
  np <- modules$np
  skimage <- modules$skimage
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


#' Function to find the centreline of each scar in mapped scar pixel data
#' @param scarData Matrix of zeros and ones, where the ones are the pixels with a scar on them
#' @param np numpy from package reticulate, or NULL to import it
#' @param skimage skimage.morphology from package reticulate, or NULL to import it
#' @param periodic Set to TRUE if data is generated with periodic boundary conditions
#'
#' @returns List with the centreline of each scar, as a data frame of pixel coordinates. The list is empty if there are
#' no scars
#' @export
find.scars.centrelines.specific.time <- function(scarData, np = NULL, skimage = NULL, periodic = TRUE){
  modules <- python.modules(np = np, skimage = skimage)
  components <- find.connected.components(mask = scarData, periodic = periodic)
  return(lapply(components, skeletonize_object_py, np = modules$np, skimage = modules$skimage))
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
  # Applies fun to the dataframe of every timestep, skipping empty timesteps. Empty timesteps stored as 0 by older
  # versions of the package become data frames with zero rows, as find.vortices.over.time and find.scars.over.time
  # give now
  apply.to.points <- function(points, fun){
    if(is.data.frame(points)){
      return(fun(points))
    }
    for(i in seq_along(points)){
      if(is.numeric(points[[i]])){
        points[[i]] <- structures.to.data.frame(list())
      } else if(!is.empty.snapshot(points[[i]])){
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
