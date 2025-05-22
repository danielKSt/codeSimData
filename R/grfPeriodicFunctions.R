#' Function to simulate a two dimensional gaussian random field with periodic boundary conditions using the SPDE approach
#'
#' @param discParamX Discretization step size in X-direction
#' @param discParamY Discretization step size in Y-direction. If none is provided, discParamX is used for both X and Y axis
#' @param theta Parameter values, (tau, kappa)
#' @param limX How far to go in x-direction
#' @param limY How far to go in y-direction
#' @param torus Mesh to use for simulation
#' @param femMatrices precision matrices, output from fmesher::fm_fem(torus)
#' @param plotResult Set to TRUE if you want a plot of the resulting
#' @param outputAsDataframe If set to true, the output is a dataframe with coordinates and values as columns
#'
#'
#' @export
simu.periodic.grf.2d <- function(discParamX, theta, limX = 1.0, limY = 1.0, discParamY = NULL, torus = NULL, femMatrices = NULL, plotResult = FALSE, outputAsDataframe = FALSE){
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
  Q <- theta[1]^2*(theta[2]^4*femMatrices$cc+2*theta[2]^2*femMatrices$g1+femMatrices$g2)
  z <- stats::rnorm(n = dim(Q)[1])
  L <- Matrix::Cholesky(A = Q, LDL = FALSE, perm = TRUE)
  y <- as.vector(Matrix::solve(L, Matrix::solve(L, z, system = 'Lt'), system = 'Pt'))
  if(plotResult){
    graphics::par(mfrow = c(1,1))
    Matrix::image(matrix(data = y, nrow = discParamX, ncol = discParamY))
  }
  if(!outputAsDataframe){
    invisible(y)
  } else {
    simResDF <- data.frame(matrix(data = NA, ncol = 3, nrow = dim(Q)[1]))
    colnames(simResDF)[3] <- "SimulationResult"
    simResDF$SimulationResult <- y
    simResDF[, 1:2] <- expand.grid(X = 1:discParamX, Y = 1:discParamY)
    simResDF$X1 <- simResDF$X1*limX/discParamX
    simResDF$X2 <- simResDF$X2*limY/discParamY
    invisible(simResDF)
  }
}

#' Function to estimate parameters for a periodic GRF
#'
#' @param observedData data to use for estimation
#' @param priorTheta Parameters for prior
#' @param spde inla.spde2 object
#' @param torus Mesh
#' @param A A-matrix
#' @param femMatrices FEM approximation matrices
#' @param reducedOutput Set to TRUE if you want only summary.hyperpar and marginal.hyperpar
#'
#' @returns INLA output from estimation of grf parameters
#'
#' @export

estimate.param.periodic.grf <- function(observedData, priorTheta = c(1, 1, 1), spde = NULL, torus = NULL, A = NULL, femMatrices = NULL, reducedOutput = FALSE){
  if(is.null(torus) && is.null(A)){
    return("Error: Not possible to make A-matrix")
  } else if(is.null(A)) {
    A <- fmesher::fm_basis(torus, loc = list(observedData$X1, observedData$X2))
  }

  if(is.null(spde)){
    if(is.null(femMatrices) && is.null(torus)){
      return("Error: Not possible to make spde object")
    } else if(is.null(femMatrices)){
      femMatrices <- fmesher::fm_fem(torus)
    }
    spde <- INLA::inla.spde2.generic(M0 = femMatrices$cc, M1 = femMatrices$g1, M2 = femMatrices$g2,
                               B0 = matrix(c(0, 1, 0), 1, 3), B1 = matrix(c(0, 0, 1), 1, 3),
                               B2 = 1, theta.mu = c(log(priorTheta[1]),log(priorTheta[2])),
                               theta.Q = diag(x = priorTheta[3], nrow = 2, ncol = 2), transform = "identity")
  }
  stk <- INLA::inla.stack(
    data = list(resp = observedData$SimulationResult),
    A = list(A, 1),
    effects = list(i = 1:spde$n.spde,
                   beta0 = rep(1, nrow(observedData))),
    tag = 'est'
  )

  res <- INLA::inla(resp ~ 0 + f(i, model = spde),
              data = INLA::inla.stack.data(stk),
              control.predictor = list(A = INLA::inla.stack.A(stk)))

  if(reducedOutput){
    invisible(list(summary.hyperpar = res$summary.hyperpar, marginals.hyperpar = res$marginals.hyperpar))
  } else {
    invisible(res)
  }
}

#' Function for plotting the theta values as wanted
#' @param res output from parameter estimation
#' @param trueParams True parameter values
#' @param kappaUpperLim upper limit on the kappa parameter for plotting range
#'
#' @export
plotting.transformed.theta <- function(res, trueParams, kappaUpperLim = NULL){
  marginalTransformedTau <- data.frame(data = res$marginals.hyperpar[[2]])
  marginalTransformedTau$data.y <- marginalTransformedTau$data.y/exp(marginalTransformedTau$data.x)
  marginalTransformedTau$data.x <- exp(marginalTransformedTau$data.x)

  marginalTransformedKappa <- data.frame(data = res$marginals.hyperpar[[3]])
  marginalTransformedKappa$data.y <- marginalTransformedKappa$data.y/exp(marginalTransformedKappa$data.x)
  marginalTransformedKappa$data.x <- exp(marginalTransformedKappa$data.x)

  if(is.null(kappaUpperLim)){kappaUpperLim <- max(c(marginalTransformedKappa$data.x), trueParams[2])}
  graphics::par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
  plot(marginalTransformedTau, type = "l",
       xlab = expression(tau), ylab = 'Posterior density')
  graphics::abline(v = exp(res$summary.hyperpar$`0.5quant`[2]), col = "red")
  graphics::abline(v = trueParams[1], col = "purple")
  plot(marginalTransformedKappa, type = "l", xlim = c(0, kappaUpperLim),
       xlab = expression(kappa), ylab = 'Posterior density')
  graphics::abline(v = exp(res$summary.hyperpar$`0.5quant`[3]), col = "red")
  graphics::abline(v = trueParams[2], col = "purple")
}

#' Function to estimate parameters for a periodic GRF, using PC prior
#'
#' @param observedData data to use for estimation
#' @param priorSigma Prior parameters for sigma
#' @param priorRange Prior parameters for the range
#' @param spde inla.spde2 object
#' @param torus Mesh
#' @param A A-matrix
#' @param femMatrices FEM approximation matrices
#' @param reducedOutput Set to TRUE if you want only summary.hyperpar and marginal.hyperpar
#'
#' @returns INLA output from estimation of grf parameters
#'
#' @export

estimate.param.periodic.grf.pcprior <- function(observedData, priorSigma = c(1, 1), priorRange = c(1,1), spde = NULL, torus = NULL, A = NULL, femMatrices = NULL, reducedOutput = FALSE){
  if(is.null(torus) && is.null(A)){
    return("Error: Not possible to make A-matrix")
  } else if(is.null(A)) {
    A <- fmesher::fm_basis(torus, loc = list(observedData$X1, observedData$X2))
  }

  if(is.null(spde)){
    if(is.null(femMatrices) && is.null(torus)){
      return("Error: Not possible to make spde object")
    } else if(is.null(femMatrices)){
      femMatrices <- fmesher::fm_fem(torus)
    }
    spde <- INLA::inla.spde2.generic(M0 = femMatrices$cc, M1 = femMatrices$g1, M2 = femMatrices$g2,
                                     B0 = cbind(-0.5*log(4*pi) - log(8)/2, 1, -1), B1 = cbind(log(8)/2, -1, 0),
                                     B2 = 1, theta.mu = c(log(1),log(1)),
                                     theta.Q = diag(x = 1, nrow = 2, ncol = 2), transform = "identity")

    lam1 <- -log(priorRange[2]) * priorRange[1]
    initial.range <- log(priorRange[1]) + 1

    lam2 <- -log(priorSigma[2]) / priorSigma[1]
    initial.sigma <- log(priorSigma[1]) - 1

    pcmatern.param <- c(lam1, lam2, 2)

    spde$f$hyper.default <-
      list(
        theta1 = list(
          prior = "pcmatern",
          param = pcmatern.param,
          initial = initial.range,
          fixed = FALSE
        ),
        theta2 = list(
          initial = initial.sigma,
          fixed = FALSE
        )
      )

  }
  stk <- INLA::inla.stack(
    data = list(resp = observedData$SimulationResult),
    A = list(A, 1),
    effects = list(i = 1:spde$n.spde,
                   beta0 = rep(1, nrow(observedData))),
    tag = 'est'
  )

  res <- INLA::inla(resp ~ 0 + f(i, model = spde),
                    data = INLA::inla.stack.data(stk),
                    control.predictor = list(A = INLA::inla.stack.A(stk)))

  if(reducedOutput){
    invisible(list(summary.hyperpar = res$summary.hyperpar, marginals.hyperpar = res$marginals.hyperpar))
  } else {
    invisible(res)
  }
}


# A simple showcase of how one might use the functions from this file
# torus <- fmesher::fm_tensor(list(fmesher::fm_mesh_1d(loc = seq(from = 0, to = 1, by = 1/64), interval = c(0,1), boundary = "cyclic"),
#                                  fmesher::fm_mesh_1d(loc = seq(from = 0, to = 1, by = 1/64), interval = c(0,1), boundary = "cyclic")))
# a <- simu.periodic.grf.2d(discParamX = 256, theta = c(0.1, 5), outputAsDataframe = TRUE, plotResult = FALSE)
# b <- estimate.param.periodic.grf(observedData = a, priorTheta = c(0.1, 5, 10), torus = torus, reducedOutput = TRUE)
# plotting.transformed.theta(res = b, trueParams = c(0.1, 5))
