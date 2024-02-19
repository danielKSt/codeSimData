#### Load data, libraries and compiled c-code  ----

## Libraries

library(tidyverse)
library("rhdf5")

## Data

setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData")
completeSurfaceDivergenceData <- h5read(file = 'Data/surface_vortices.h5', name = 'sdiv')

## C-code

dyn.load("surfaceExploration.so")

#### Helper functions for exploration ----

distance <- function(x,y){
  return(sqrt(min((x[1]-y[1])^2, (x[1]+1-y[1])^2, (y[1]+1-x[1])^2)+min((x[2]-y[2])^2, (x[2]+1-y[2])^2, (y[2]+1-x[2])^2)))
}

make.grid.of.distances <- function(locs, metric.function){
  N = dim(locs)[1]
  distances <- matrix(nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in i:N) {
      distances[i,j] <- metric.function(locs[i,], locs[j,])
      distances[j,i] <- distances[i,j]
    }
  }
  return(distances)
}

estimate.variogram <- function(measurementLocations, measurements, bins, metric.function){
  distances <- make.gird.of.distances(measurementLocations, metric.function)
  estimatedVariogram <- c(1:length(bins))
  for (k in c(1:length(bins))) {
    
  }
}


#### Exploration of surface divergence ----
