
# Testing for data preparation

testVortex <- matrix(data = c(0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0), nrow = 5)
testingPixels <- list(c(2, 2), c(2, 3), c(3, 1), c(3, 2), c(3, 3), c(4, 2))

explore.vortex <- function(vortex.data, x_start, y_start){
  remaining <- list(c(x_start, y_start))
  remainingCount <- 1
  vortex.data[x_start, y_start] <- 2
  pixels <- data.frame(x = x_start, y = y_start)
  pixelsCount <- 1
  while(remainingCount > 0){
    x <- remaining[[remainingCount]][1]
    y <- remaining[[remainingCount]][2]
    remaining[[remainingCount]] <- NULL
    remainingCount <- remainingCount - 1
    kj_list <- list(c(1, -1), c(1, 0), c(-1, 1), c(0, 1), c(1, 1))
    for (kj in kj_list) {
      k <- kj[1]
      j <- kj[2]
      if(x + k > 0 && x + k < dim(vortex.data)[1] && y + j > 0 && y + j < dim(vortex.data)[2] && vortex.data[x+k, y+j] == 1){
        remainingCount <- remainingCount + 1
        remaining[[remainingCount]] <- c(x+k, y+j)
        vortex.data[x+k, y+j] <- 2
        pixels[nrow(pixels)+1,] <- c(x+k, y+j)
        pixelsCount <- pixelsCount + 1
      }
    }
  }
  return(pixels)
}


vortex.data <- testVortex
x_start <- 2
y_start <- 2

resultingPixels <- explore.vortex(vortex.data = vortex.data, x_start = x_start, y_start = y_start)

calculate.vortex.properties <- function(pixels){
  centreCoordinates <- c(round(mean(pixels$x)), round(mean(pixels$y)))
  noPixels <- dim(pixels)[1]
  xStretch <- max(pixels$x)-min(pixels$x)
  yStretch <- max(pixels$y)-min(pixels$y)
  return(c(round(mean(pixels$x)), round(mean(pixels$y)), dim(pixels)[1], max(pixels$x)-min(pixels$x), max(pixels$y)-min(pixels$y)))
}

