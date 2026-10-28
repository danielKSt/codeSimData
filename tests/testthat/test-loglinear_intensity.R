# Number of points in each pixel, with the pixel convention of get_Z_at_coord
pixel.counts <- function(Z, structures, l_x, l_y){
  counts <- array(data = 0L, dim = dim(Z))
  for(t in seq_along(structures)){
    for(i in seq_len(nrow(structures[[t]]))){
      xInd <- min(floor(structures[[t]]$x[i]*dim(Z)[2]/l_x) + 1, dim(Z)[2])
      yInd <- min(floor(structures[[t]]$y[i]*dim(Z)[3]/l_y) + 1, dim(Z)[3])
      counts[t, xInd, yInd] <- counts[t, xInd, yInd] + 1L
    }
  }
  return(counts)
}

# With a covariate that is constant within each pixel, the Poisson process likelihood is the likelihood of a Poisson
# GLM of the pixel counts with offset log(pixelArea), except for the terms n*log(pixelArea) and sum(log(counts!))
glm.logLik.as.process <- function(glmFit, counts, pixelArea){
  return(as.numeric(stats::logLik(glmFit)) - sum(counts)*log(pixelArea) + sum(lfactorial(counts)))
}

test_that("The model on log scale gives the estimates of a Poisson GLM of the pixel counts", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  field <- setup$field[tidspunkt, , , drop = FALSE]
  pixelArea <- setup$l_x*setup$l_y/(dim(field)[2]*dim(field)[3])

  # The field itself has negative values, and its positive part has zeros
  for(Z in list(field, pmax(field, 0))){
    for(structures in points){
      fit <- fit.given.covar.loglinear(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                       returnParams = TRUE)
      expect_equal(fit.given.covar.loglinear(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                             returnParams = FALSE), fit$logLik)

      counts <- pixel.counts(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y)
      glmFit <- stats::glm(as.vector(counts) ~ as.vector(Z), family = stats::poisson(),
                           offset = rep(log(pixelArea), length(Z)), control = stats::glm.control(epsilon = 1e-12))
      # optimize finds beta_1 to about 1e-4 relative, as for fit.given.covar, while the log likelihood is flat there
      expect_equal(unname(stats::coef(glmFit)), c(fit$beta0, fit$beta1), tolerance = 1e-4)
      expect_equal(fit$logLik, glm.logLik.as.process(glmFit, counts, pixelArea), tolerance = 1e-8)
    }
  }
})

test_that("The power model is the model on log scale with log(Z)", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- exp(setup$field[tidspunkt, , , drop = FALSE])
  pixelArea <- setup$l_x*setup$l_y/(dim(Z)[2]*dim(Z)[3])

  for(structures in points){
    power <- fit.given.covar(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y, returnParams = TRUE)
    logScale <- fit.given.covar.loglinear(Z = log(Z), structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                          returnParams = TRUE)
    expect_equal(logScale$logLik, power$logLik)
    expect_equal(logScale$beta1, power$b, tolerance = 1e-5)
    expect_equal(exp(logScale$beta0), power$a, tolerance = 1e-5)

    counts <- pixel.counts(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y)
    glmFit <- stats::glm(as.vector(counts) ~ log(as.vector(Z)), family = stats::poisson(),
                         offset = rep(log(pixelArea), length(Z)), control = stats::glm.control(epsilon = 1e-12))
    expect_equal(power$logLik, glm.logLik.as.process(glmFit, counts, pixelArea), tolerance = 1e-8)
  }
})

test_that("The fit does not depend on the units or the origin of the covariate", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- pmax(setup$field[tidspunkt, , , drop = FALSE], 0)
  fit <- fit.given.covar.loglinear(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                   returnParams = TRUE)

  for(scale in c(1e-8, 1e8)){
    scaled <- fit.given.covar.loglinear(Z = scale*Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                        returnParams = TRUE)
    expect_equal(scaled$logLik, fit$logLik)
    expect_equal(scale*scaled$beta1, fit$beta1, tolerance = 1e-5)
  }
  shifted <- fit.given.covar.loglinear(Z = Z + 100, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                       returnParams = TRUE)
  expect_equal(shifted$logLik, fit$logLik)
  expect_equal(shifted$beta1, fit$beta1, tolerance = 1e-5)
})

test_that("A constant covariate gives the constant intensity", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- array(data = 0.3, dim = c(length(tidspunkt), dim(setup$field)[2:3]))
  n <- sum(vapply(points$scars, nrow, integer(1)))
  totalArea <- length(tidspunkt)*setup$l_x*setup$l_y

  fit <- fit.given.covar.loglinear(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                   returnParams = TRUE)
  expect_equal(fit$beta1, 0)
  expect_equal(fit$logLik, n*log(n/totalArea) - n)
})

test_that("Points only where the covariate is largest give the limit of an infinite beta_1", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  # 1 on the pixels with scars and 0 elsewhere
  Z <- pmin(pixel.counts(Z = array(data = 0, dim = c(length(tidspunkt), dim(setup$field)[2:3])),
                         structures = points$scars, l_x = setup$l_x, l_y = setup$l_y), 1)
  n <- sum(vapply(points$scars, nrow, integer(1)))
  pixelArea <- setup$l_x*setup$l_y/(dim(Z)[2]*dim(Z)[3])

  # The supremum of the likelihood is a constant intensity on the pixels with scars and zero elsewhere
  fit <- fit.given.covar.loglinear(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                   returnParams = TRUE)
  expect_gt(fit$beta1*stats::sd(as.vector(Z)), 5)
  expect_equal(fit$logLik, n*log(n/(sum(Z)*pixelArea)) - n)
})

test_that("Covariates that are not finite give an error", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- setup$field[tidspunkt, , , drop = FALSE]
  for(value in c(NA, NaN, Inf, -Inf)){
    Z[1, 1, 1] <- value
    expect_error(fit.given.covar.loglinear(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                           returnParams = FALSE), "must be finite")
  }
})
