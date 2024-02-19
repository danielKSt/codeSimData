#### Set wd and load libraries and code ----

library(testthat)

## Working directory
setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData")

## C-code

dyn.unload("surfaceExploration.so")
dyn.load("surfaceExploration.so")

#### Tests ----

test.TorusDistance <- function(tol){
  d_1 <- 1
  x_1 <- c(0.4)
  y_1 <- c(0.3)
  distance <- .C("distanceOnTorus", as.integer(d_1), as.double(x_1), as.double(y_1), result = as.double(0))
  if((distance[['result']]-0.1) < tol){
    cat("Success for first test")
  }
}

test.TorusDistance(1e-10)

