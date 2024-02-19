#### Load data  ----

library("rhdf5")
setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Simulated data")
completeSurfaceVortexData <- h5read(file = 'surface_vortices.h5', name = 'surf')
completeSurfaceDivergenceData <- h5read(file = 'surface_vortices.h5', name = 'sdiv')

library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(reshape)

endTimestep <- dim(completeSurfaceVortexData)[1]


#### Functions for data preparation into point process format ----

explore.vortex <- function(vortexPixels, boxDimensions){
  if(any(vortexPixels$x==1) && any(vortexPixels$x==boxDimensions[1])){
    vortexPixels$x <- vortexPixels$x*(vortexPixels < floor(boxDimensions[1]/2)) - (-1*vortexPixels$x %% boxDimensions[1])*(vortexPixels$x >= floor(boxDimensions[1]/2))
    xMean <- mean(vortexPixels$x)%% boxDimensions[1]
  }
  if(any(vortexPixels$y==1) && any(vortexPixels$y==boxDimensions[2])){
    vortexPixels$y <- vortexPixels$y*(vortexPixels < floor(boxDimensions[2]/2)) - (-1*vortexPixels$y %% boxDimensions[2])*(vortexPixels$y >= floor(boxDimensions[2]/2)) 
    yMean <- mean(vortexPixels$y)%% boxDimensions[2]
  }
  return(c(xMean/boxDimensions[2], yMean/boxDimensions[1], dim(vortexPixels$x), max(vortexPixels$x)-min(vortexPixels$x), max(vortexPixels$y)-min(vortexPixels$y)))
}

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
        vortexSummary <- c(mean(pixels$x)/256.0, mean(pixels$y)/256.0, dim(pixels)[1], max(pixels$x)-min(pixels$x), max(pixels$y)-min(pixels$y))
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

# We make a list with the point process data
vortices <- vector(mode = "list", length = endTimestep)

for(i in 1:endTimestep) {
  vortices[[i]] <- find.vortices.specific.time(completeSurfaceVortexData[i, 1:256, 1:256])
}

numberOfVortices <- c(1:endTimestep)
pixelsWithVortex <- c(1:endTimestep)
betaSquared <- c(1:endTimestep)
for (i in c(1:endTimestep)) {
  if(!is.numeric(vortices[[i]])){
    numberOfVortices[i] <- dim(vortices[[i]])[1]
  } else {
    numberOfVortices[i] <- 0
  }
  pixelsWithVortex[i] <- sum(completeSurfaceVortexData[i, 1:256, 1:256])
  betaSquared[i] <- sum(completeSurfaceDivergenceData[i, 1:256, 1:256]^2)
}
numberOfVortices_rescaled <- (100/max(numberOfVortices))*numberOfVortices
pixelsWithVortex <- (100/max(pixelsWithVortex))*pixelsWithVortex
betaSquared_rescaled <- (100/max(betaSquared))*betaSquared

numberOfVorticesAndBetaSquared <- data.frame(t = c(1:endTimestep), noVortices = numberOfVortices_rescaled, betaSq = betaSquared_rescaled)
numberOfVorticesAndBetaSquared <- melt(numberOfVorticesAndBetaSquared, id.vars = "t")
a <- ggplot(data = numberOfVorticesAndBetaSquared, aes(x = t, y = value, color = variable))+
  geom_line()

a

#### We do some simple plotting at some specific snapshots in time ----

timeOfSnapshot <- c(100, 2000, 5500, 8500)
#timeOfSnapshot <- c(2000)

makePlotVorticeBetaAtTime <- function(timeIndex, betaData, vortexData){
  vortexPlot <- ggplot(data = vortexData[[timeIndex]], mapping = aes(x = x, y = y, size = noPixels)) +
    geom_point(alpha = 0.7) +
    theme(legend.position = "bottom")+
    ggtitle(sprintf("Vortices at time %i", timeIndex))
  
  data <- expand.grid(X = c(1:256)/256.0, Y = c(1:256)/256.0)
  data$SurfaceDivergence <- array(betaData[timeIndex, 1:256, 1:256])
  betaPlot <- ggplot(data, mapping = aes(x = X, y = Y, fill = SurfaceDivergence))+
    geom_tile() +
    scale_size(range = c(1.0, 10), name = "Area of vortex (Pixels)") +
    theme(legend.position = "bottom")+
    ggtitle(sprintf("Surface divergence at time %i", timeIndex))
  return(list("betaPlot" = betaPlot, "vortexPlot" = vortexPlot))
}

vortexAndBetaPlot <- timeOfSnapshot |> lapply(FUN = makePlotVorticeBetaAtTime, betaData = completeSurfaceDivergenceData, vortexData = vortices)

library("gridExtra")
png("Presentation_december23/vortexPlots1.png", res = 400, width = 10000, height = 6179)
grid.arrange(vortexAndBetaPlot[[1]]$vortexPlot, vortexAndBetaPlot[[1]]$betaPlot, vortexAndBetaPlot[[2]]$vortexPlot, 
             vortexAndBetaPlot[[2]]$betaPlot, ncol = 2, nrow = 2)
dev.off()
png("Presentation_december23/vortexPlots2.png", res = 400, width = 10000, height = 6179)
grid.arrange(vortexAndBetaPlot[[3]]$vortexPlot, vortexAndBetaPlot[[3]]$betaPlot, vortexAndBetaPlot[[4]]$vortexPlot, 
             vortexAndBetaPlot[[4]]$betaPlot, ncol = 2, nrow = 2)
dev.off()


#### We do some exploratory data analysis of the point process data ----

# We look at intensity
stationaryIntensity <- mean(numberOfVortices)
poissonQuantiles <- qpois(p = c(0.025, 0.975), lambda = stationaryIntensity)
sum((numberOfVortices < poissonQuantiles[1]))/endTimestep
sum((numberOfVortices > poissonQuantiles[2]))/endTimestep

calcAreaOfBallInSquare <- function(centre, radius){
  m = 10000
  pointsInBall <- centre + t(circleUnif(n = m, r = radius))
  pointsInSquareAndBall <- (pointsInBall[1, 1:m] > 0.0)*(pointsInBall[1, 1:m] < 1.0) * (pointsInBall[2, 1:m] > 0.0)*(pointsInBall[2, 1:m] < 1.0)
  return((sum(pointsInSquareAndBall)/m)*pi*radius^2)
}

calcLambdaFromRadiusAndLocation <- function(time, vortices, radius, location){
  noVorticesInBall <- 0
  for(i in c(1:length(vortices[[time]]))){
    if(((vortices[[time]]$x[i]-location[1])^2+(vortices[[time]]$x[i]-location[1])^2)< radius^2){
      noVorticesInBall <- noVorticesInBall + 1
    }
  }
  areaOfBallInSquare <- calcAreaOfBallInSquare(centre = location, radius = radius)
  return(noVorticesInBall/areaOfBall)
}

# We look at the intensity for four snapshots
timeOfSnapshot <- c(100, 2000, 5500, 8500)
radii <- c(0.05, 0.15, 0.25, 0.5, 0.8, 1.0)
locs <- matrix(data = c(0.25, 0.25, 0.25, 0.75, 0.75, 0.25, 0.75, 0.75), ncol = 2, byrow = TRUE)
lambda <- vector(mode = "list", length = length(timeOfSnapshot))
for (i in c(1:length(timeOfSnapshot))) {
  lambda[[i]] <- matrix(ncol = dim(locs)[1], nrow = length(radii))
  for(j in c(1:length(radii))){
    for(k in c(1:dim(locs)[1])){
      lambda[[i]][j, k] <- calcLambdaFromRadiusAndLocation(time = timeOfSnapshot[i], vortices = vortices, radius = radii[j], location = locs[k])
    }
  }
}
lambda[[1]]
lambda[[2]]
lambda[[3]]
lambda[[4]]







data <- expand.grid(X = 1:256, Y = 1:256)
data$vorticesOnPixel <- array(completeSurfaceVortexData[1500, 1:256, 1:256])
ggplot(data, mapping = aes(x = X, y = Y, fill = vorticesOnPixel))+
  geom_tile()




