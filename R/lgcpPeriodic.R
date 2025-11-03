#' Function to simulate inhomogeneous poisson point process
#'
#' @param lambda Intensity of the field
#' @param discParamX Discretization step size in X-direction
#' @param discParamY Discretization step size in Y-direction. If none is provided, discParamX is used for both X and Y axis
#' @param limX How far to go in x-direction
#' @param limY How far to go in y-direction
#'
#' @importFrom dplyr filter select
#' @importFrom rlang .data
#' @export
simu.pp.from.intensity <- function(lambda, discParamX, discParamY = NULL, limX = 1, limY = 1){
  if(is.null(discParamY)){
    discParamY <- discParamX
  }
  nPointsUniform <- stats::rpois(n = 1, lambda = max(lambda))
  initialPPP <- data.frame(matrix(data = stats::runif(n = 4*nPointsUniform), nrow = nPointsUniform, ncol = 4))
  initialPPP$X1 <- initialPPP$X1*limX
  initialPPP$X2 <- initialPPP$X2*limY
  for (i in 1:nPointsUniform) {

    xInd <- round(initialPPP$X1[i]*(discParamX/limX))+(discParamX/limX)*(round(initialPPP$X1[i]*(discParamX/limX))==0)
    yInd <- round(initialPPP$X2[i]*(discParamY/limY))+(discParamY/limY)*(round(initialPPP$X2[i]*(discParamY/limY))==0)
    initialPPP$X4[i] <- lambda[xInd,yInd]/max(lambda)
  }
  PPP <- initialPPP |>
    dplyr::filter(.data$X4 > .data$X3) |>
    dplyr::select("X1", "X2")
  colnames(PPP) <- c("x", "y")
  return(PPP)
}



#' Function to simulate a log-Gaussian Cox process with periodic boundary conditions
#'
#' @param discParamX Discretization step size in X-direction
#' @param discParamY Discretization step size in Y-direction. If none is provided, discParamX is used for both X and Y axis
#' @param theta Parameter values, (tau, kappa)
#' @param beta lambda = exp(beta1+beta2*grf)
#' @param limX How far to go in x-direction
#' @param limY How far to go in y-direction
#' @param torus Mesh to use for simulation
#' @param femMatrices precision matrices, output from fmesher::fm_fem(torus)
#' @param plotLambda Set to TRUE if you want a plot of the simulated intensity
#' @param plotVortices Set to TRUE if you want a plot of the simulated vortices
#' @param includeLambdaInOutput Set to false if you only want the vortices as output
#'
#' @export
simu.lgcp.pp <- function(discParamX, beta, theta, limX = 1.0, limY = 1.0, discParamY = NULL, torus = NULL, femMatrices = NULL, plotLambda = FALSE, plotVortices = FALSE, includeLambdaInOutput = TRUE){
  if(is.null(discParamY)){
    discParamY <- discParamX
  }
  if(is.null(femMatrices) && is.null(torus)){
    torus <- fmesher::fm_tensor(list(fmesher::fm_mesh_1d(loc = seq(from = 0, to = limX, by = limX/discParamX), interval = c(0,limX), boundary = "cyclic"),
                                     fmesher::fm_mesh_1d(loc = seq(from = 0, to = limY, by = limY/discParamY), interval = c(0,limY), boundary = "cyclic")))
    femMatrices <- fmesher::fm_fem(torus)
  } else if(is.null(femMatrices)) {
    femMatrices <- fmesher::fm_fem(torus)
  }
  grf <- matrix(data = simu.periodic.grf.2d(discParamX = discParamX, torus = torus, femMatrices = femMatrices, theta = theta), ncol = discParamX, nrow = discParamY)
  lambda <- exp(beta[2]*grf+beta[1])
  if(plotLambda){Matrix::image(lambda)}
  vortices <- simu.pp.from.intensity(lambda = lambda, discParamX = discParamX, discParamY = discParamY, limX = limX, limY = limY)
  if(plotVortices){plot(x = vortices[ , 1], y = vortices[ , 2], xlim = c(0, 1), ylim = c(0,1))}
  if(includeLambdaInOutput){
    return(list(vortices = vortices, lambda = lambda))
  } else {
    return(vortices)
  }
}
