#' Functions used for manipulation of the surface divergence


#' Function for estimating the norm of the gradient field of the surface divergence
#'
#' @param sDiv surface divergence
#' @param xdim Dimension in the x direction
#' @param ydim Dimension in the y direction
#' @param hx Step size in x-direction
#' @param hy Step size in y-direction
#'
#' @returns Matrix with the norm of the gradient of the surface divergence
#'
#' @export
estimate.gradient.norm.sdiv <- function(sDiv, xdim, ydim, hx, hy){
  res <- matrix(data = NA, nrow = xdim, ncol = ydim)
  # Corner pixels
  res[1,1] <- sqrt(((sDiv[2,1]-sDiv[xdim,1])/(2*hx))^2+((sDiv[1, 2]-sDiv[1,ydim])/(2*hy))^2)
  res[xdim,1] <- sqrt(((sDiv[1,1]-sDiv[xdim-1,1])/(2*hx))^2+((sDiv[xdim, 2]-sDiv[xdim,ydim])/(2*hy))^2)
  res[1,ydim] <- sqrt(((sDiv[2,ydim]-sDiv[xdim,ydim])/(2*hx))^2+((sDiv[1, 1]-sDiv[1,ydim-1])/(2*hy))^2)
  res[xdim,ydim] <- sqrt(((sDiv[1,ydim]-sDiv[xdim-1,ydim])/(2*hx))^2+((sDiv[xdim, 1]-sDiv[xdim,ydim-1])/(2*hy))^2)

  # Edge rows
  for (k in 2:(xdim-1)) {
    res[k, 1] <- sqrt(((sDiv[k+1,1]-sDiv[k-1,1])/(2*hx))^2+((sDiv[k,2]-sDiv[k,ydim])/(2*hy))^2)
    res[k, ydim] <- sqrt(((sDiv[k+1,ydim]-sDiv[k-1,ydim])/(2*hx))^2+((sDiv[k,1]-sDiv[k,ydim-1])/(2*hy))^2)
  }

  for (k in 2:(ydim-1)) {
    res[1, k] <- sqrt(((sDiv[2, k]-sDiv[xdim,k])/(2*hx))^2+((sDiv[1,k+1]-sDiv[1,k-1])/(2*hy))^2)
    res[xdim, k] <- sqrt(((sDiv[1, k]-sDiv[xdim-1,k])/(2*hx))^2+((sDiv[xdim,k+1]-sDiv[xdim,k-1])/(2*hy))^2)
  }

  # Central pixels
  for (i in 2:(xdim-1)) {
    for (j in 2:(ydim-1)) {
      res[i,j] <- sqrt(((sDiv[i+1,j]-sDiv[i-1,j])/(2*hx))^2+
                         ((sDiv[i,j+1]-sDiv[i,j-1])/(2*hy))^2)
    }
  }
  return(res)
}


#' Function for smoothing out data for some field on the surface
#'
#' @param radius Radius man ynskjer å vindauga for
#'
#' @returns Ei nx2 matrise med indexar for koordinatar innanfor ein sirkel med den gitte radiusen
#'
#' @export
finn.vindauga.sirkel <- function(radius){
  if(radius < 1){
    return(matrix(data = c(0, 0), nrow = 1, ncol = 2))
  }
  top <- round(radius)
  interval <- c(-top:top)
  vindauga <- matrix(data = NA, nrow = length(interval)^2, ncol = 2)
  accepted <- 1
  for (i in interval) {
    for (j in interval) {
      if(i*i+j*j<=radius^2){
        vindauga[accepted, ] <- c(i, j)
        accepted <- accepted + 1
      }
    }
  }
  vindauga <- vindauga[!rowSums(is.na(vindauga)),]
  return(vindauga)
}


#' Function for finding shell of a given window
#'
#' @param vindauga Vindauga me ynskjer å finne skallet til
#'
#' @returns Ei nx2 matrise med indexar for koordinatar innanfor ein sirkel med den gitte radiusen
#'
#' @export
finn.skall.av.vindauga <- function(vindauga){
  if(dim(vindauga)[1]==5){
    return(matrix(data = c(-1L, 0L, 0L, 1L, 1L, 0L, 0L, -1L), ncol = 2, byrow = TRUE))
  }
  # Initialiserer skallet
  skallDim <- 0
  skall <- matrix(data = 0, nrow = dim(vindauga)[1], ncol = dim(vindauga)[2])

  # Leggjer til pikslar frå lengst til venstre i øvre halvplan
  x <- min(vindauga[, 1])
  a <- which(vindauga[, 1] == x)
  b <- which(vindauga[a, 2] >= 0)
  indLeft <- a[b]
  for (i in indLeft) {
    skallDim <- skallDim + 1
    skall[skallDim, ] <- vindauga[i, ]
  }
  ymax <- max(vindauga[a, 2])

  for(x in (min(vindauga[, 1])+1):(-1)){
    a <- which(vindauga[, 1] == x)
    b <- which(vindauga[a, 2] > ymax)
    if(length(b) == 0){
      b <- which(vindauga[a, 2] >= ymax)
    }
    indLeft <- a[b]
    for(i in indLeft){
      skallDim <- skallDim + 1
      skall[skallDim, ] <- vindauga[i, ]
    }
    ymax <- max(vindauga[a, 2])
  }
  kvadrantDim <- skallDim
  for(i in 1:kvadrantDim){
    skall[i+kvadrantDim, ] <- c(skall[i, 2], -skall[i, 1])
    skall[i+2*kvadrantDim, ] <- c(-skall[i, 1], -skall[i, 2])
    skall[i+3*kvadrantDim, ] <- c(-skall[i, 2], skall[i, 1])
    skallDim <- skallDim+3
  }
  return(matrix(as.integer(skall[1:skallDim, ]), ncol = 2))
}

#' Funciton for smoothing out a single datapoint
#' @param sField Matrix of values on the surface
#' @param vindauga indices of sliding window relative to the pixle we smooth around
#' @param glatter matrix of weights for how to smooth the pixels
#' @param i First index
#' @param j Second index
#' @param boundaryCondition "periodic" if we have periodic data, "physical" if not
#'
#' @returns A smoothed value
single.pixle.smoothing<- function(sField, vindauga, glatter, i, j, boundaryCondition){
  xdim <- dim(sField)[1]
  ydim <- dim(sField)[2]
  if(boundaryCondition == "periodic"){
    vindaugaJustert <- t(t(vindauga) + c(i,j)) - 1
    vindaugaJustert[,1] <- vindaugaJustert[,1]%%xdim + 1
    vindaugaJustert[,2] <- vindaugaJustert[,2]%%ydim + 1
    vindaugaJustert <- lapply(c(1:dim(vindaugaJustert)[1]), function(x){return(vindaugaJustert[x,])})
  }
  else if(boundaryCondition == "physical"){
    vindaugaJustert <- t(t(vindauga) + c(i,j))
    utfor1 <- which(vindaugaJustert[ , 1] < 1)
    utfor2 <- which(vindaugaJustert[ , 2] < 1)
    utfor3 <- which(vindaugaJustert[ , 1] > xdim)
    utfor4 <- which(vindaugaJustert[ , 2] > ydim)
    utfor <- unique(c(utfor1, utfor2, utfor3, utfor4))
    if(!length(utfor)==0){vindaugaJustert <- vindaugaJustert[-utfor, ]}
    vindaugaJustert <- lapply(c(1:dim(vindaugaJustert)[1]), function(x){return(vindaugaJustert[x,])})
  }
  smoother.input <- sapply(X = vindaugaJustert,
                           function(index, sField){return(sField[index[1], index[2]])},
                           sField = sField)

  return(glatter(smoother.input))
}


#' Function for smoothing out data for some field on the surface
#'
#' @param sField Matrix of values on the surface
#' @param vindauga indices of sliding window relative to the pixle we smooth around
#' @param glatter matrix of weights for how to smooth the pixels
#' @param boundaryCondition "periodic" if we have periodic data, "physical" if not
#'
#' @returns A smoothed version of sField
#'
#' @export
smooth.surface.field <- function(sField, vindauga, glatter, boundaryCondition){
  xdim <- dim(sField)[1]
  ydim <- dim(sField)[2]
  res <- matrix(data = NA, nrow = xdim, ncol = ydim)
  indices <- as.matrix(expand.grid(X = 1:xdim, Y = 1:ydim))
  res <- apply(indices, MARGIN = 1,
               FUN = function(x){
                 return(
                   single.pixle.smoothing(sField = sField,
                                          vindauga = vindauga,
                                          glatter = glatter,
                                          i = x[1],
                                          j = x[2],
                                          boundaryCondition = boundaryCondition)
                   )
                 })
  return(matrix(data = res, nrow = xdim))
}


#' Function for calculating the local variance with periodic boundary conditions
#'
#' @param sField Matrix of values on the surface
#' @param vindauga indices of sliding window relative to the pixle we smooth around
#' @param skall indices of shell of sliding window relative to the pixle we smooth around
#'
#' @returns A smoothed version of sField
#'
#' @useDynLib codeSimData
#' @importFrom Rcpp evalCpp
#'
#' @export
locvar.transformation.periodic <- function(sField, vindauga, skall){
  if(dim(vindauga)[1]!=2){vindauga <- t(vindauga)}
  if(dim(skall)[1]!=2){skall <- t(skall)}
  res <- calculateLocVar_periodic(inWin = vindauga, inShell = skall, insDiv = sField,
                                  inDims = c(dim(sField)[1], dim(sField)[2], dim(vindauga)[2], dim(skall)[2]))
  return(matrix(data = res, nrow = dim(sField)[1]))
}


# sField <- matrix(data = c(1,2,3,4,5,6,7,
#                           2,3,4,5,6,7,8,
#                           3,4,5,6,7,8,9,
#                           4,5,6,7,8,9,10,
#                           5,6,7,8,9,10,11), nrow = 7)
#
# vindauga <- finn.vindauga.sirkel(radius = 2)
# glatter <- mean
# boundaryCondition <- "physical"
#
# a <- smooth.surface.field(sField = sField, vindauga = vindauga, glatter = mean, boundaryCondition = "periodic")
# #View(a)
# image(sField)
# image(a)

#library(Rcpp)
#compileAttributes()
# sField <- matrix(data = 0, nrow = 7, ncol = 5)
# for(y in 1:dim(sField)[2]){
#   for(x in 1:dim(sField)[1]){
#     sField[x,y] <- x+y
#   }
# }
#
# vindauga <- finn.vindauga.sirkel(radius = 2)
# skall <- finn.skall.av.vindauga(vindauga = vindauga)
# res <- calculateLocVar_periodic(inWin = t(vindauga), inShell = t(skall), insDiv = sField,
#                                 inDims = c(7, 5, dim(vindauga)[1], dim(skall)[1]))
