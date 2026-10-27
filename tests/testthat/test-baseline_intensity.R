# Covariate with no zeros (local mean square) and with zeros (the positive part of the field) at the snapshots tidspunkt
make.baseline.test.covariates <- function(setup, tidspunkt){
  field <- setup$field[tidspunkt, , , drop = FALSE]
  vindauga <- t(finn.vindauga.sirkel(radius = 2))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
  positive <- array(data = NA_real_, dim = dim(field))
  for(i in seq_along(tidspunkt)){
    positive[i, , ] <- calculateLocMeanSq_periodic(inWin = vindauga, inShell = skall, insDiv = field[i, , ],
                                                   inDims = c(dim(field)[2], dim(field)[3], dim(vindauga)[2],
                                                              dim(skall)[2]))
  }
  return(list(positive = positive, withZeros = pmax(field, 0)))
}

baseline.areas <- function(Z, setup){
  pixelArea <- setup$l_x*setup$l_y/(dim(Z)[2]*dim(Z)[3])
  return(list(pixelArea = pixelArea, totalArea = pixelArea*length(Z)))
}

test_that("Without baseline the complete likelihood is the one of the model without baseline", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$positive
  areas <- baseline.areas(Z = Z, setup = setup)
  inData <- get_sumLogZ_and_n(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y)
  logZAtPoints <- log(get_Z_at_points(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y))
  expect_length(logZAtPoints, inData$n)
  expect_equal(sum(logZAtPoints), inData$sumLogZ)

  for(b in c(0, 0.7, 3)){
    expect_equal(calc_likelihood_baseline_complete(nu0 = 0, a = 2.5, b = b, logZ = log(Z), logZAtPoints = logZAtPoints,
                                                   pixelArea = areas$pixelArea, totalArea = areas$totalArea),
                 calc_likelihood_poisson_complete(logA = log(2.5), b = b, logZ = log(Z), sumLogZ = inData$sumLogZ,
                                                  n = inData$n, pixelArea = areas$pixelArea))
  }
})

test_that("The estimates maximise the complete likelihood", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  covariates <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)

  for(Z in covariates){
    areas <- baseline.areas(Z = Z, setup = setup)
    for(structures in points){
      fit <- fit.given.covar.baseline(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                      returnParams = TRUE)
      expect_equal(fit.given.covar.baseline(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                            returnParams = FALSE), fit$logLik)
      expect_gte(fit$b, 0)
      expect_true(fit$w >= 0 && fit$w <= 1)
      expect_equal(fit$lambda, fit$nu0/fit$a)

      complete <- function(par){
        calc_likelihood_baseline_complete(nu0 = par[1], a = par[2], b = par[3], logZ = log(Z[Z > 0]),
                                          logZAtPoints = log(get_Z_at_points(Z, structures, setup$l_x, setup$l_y)),
                                          pixelArea = areas$pixelArea, totalArea = areas$totalArea)
      }
      expect_equal(complete(c(fit$nu0, fit$a, fit$b)), fit$logLik, tolerance = 1e-8)
      # A direct maximisation of the complete likelihood, from several starting points, is not better. nu0 = 0 gives
      # -Inf when there are points on zero pixels, so nu0 and a are kept just above 0
      for(start in list(c(fit$nu0, fit$a, fit$b) + 0.1, c(1, 1, 0.5), c(0.1, 2, 2))){
        opt <- stats::optim(par = start, fn = complete, method = "L-BFGS-B", lower = c(1e-10, 1e-10, 0),
                            control = list(fnscale = -1, factr = 1e3))
        expect_lte(opt$value, fit$logLik + 1e-4)
      }
    }
  }
})

test_that("The model with baseline is at least as good as the model without, and the same when w = 0", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$positive

  for(structures in points){
    fit <- fit.given.covar.baseline(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                    returnParams = TRUE)
    fitPower <- fit.given.covar(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y, returnParams = TRUE)
    if(fitPower$b >= 0){
      expect_gte(fit$logLik, fitPower$logLik - 1e-8)
    }
  }

  # Three pixels with Z = 1 and one with Z = 10, with all the points on the last one. The derivative of the profile
  # likelihood is then negative at w = 0, so the estimate is exactly w = 0
  logZ <- log(c(1, 1, 1, 10))
  logZAtPoints <- rep(log(10), 5)
  fitW <- find_w_mle(b = 1, logZ = logZ, logZAtPoints = logZAtPoints, pixelArea = 0.5, totalArea = 2)
  expect_identical(fitW$w, 0)
  expect_equal(fitW$logLik, 5*log(5) - 5 + calc_likelihood_poisson_b(b = 1, logZ = logZ, sumLogZ = sum(logZAtPoints),
                                                                     n = 5, pixelArea = 0.5))
})

test_that("Zero covariates get the baseline intensity only, also for b = 0", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$withZeros
  areas <- baseline.areas(Z = Z, setup = setup)
  zAtPoints <- get_Z_at_points(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y)
  expect_true(any(zAtPoints == 0))

  # With b = 0 the intensity is nu0 on the zero pixels and nu0 + a on the others
  nu0 <- 0.3
  a <- 1.2
  expected <- sum(log(ifelse(zAtPoints > 0, nu0 + a, nu0))) - nu0*areas$totalArea - a*areas$pixelArea*sum(Z > 0)
  expect_equal(calc_likelihood_baseline_complete(nu0 = nu0, a = a, b = 0, logZ = log(Z[Z > 0]),
                                                 logZAtPoints = log(zAtPoints), pixelArea = areas$pixelArea,
                                                 totalArea = areas$totalArea), expected)

  fit <- fit.given.covar.baseline(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                  returnParams = TRUE)
  expect_true(is.finite(fit$logLik))
  expect_gt(fit$w, 0)
})

test_that("Points only on zero pixels give a constant intensity", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$withZeros
  areas <- baseline.areas(Z = Z, setup = setup)
  # One point at the centre of the first zero pixel of each snapshot
  structures <- lapply(seq_along(tidspunkt), function(t){
    ind <- which(Z[t, , ] == 0, arr.ind = TRUE)[1, ]
    data.frame(x = (ind[1] - 0.5)*setup$l_x/dim(Z)[2], y = (ind[2] - 0.5)*setup$l_y/dim(Z)[3])
  })
  fit <- fit.given.covar.baseline(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                  returnParams = TRUE)
  n <- length(tidspunkt)
  expect_identical(fit$w, 1)
  expect_identical(fit$a, 0)
  expect_equal(fit$nu0, n/areas$totalArea)
  expect_equal(fit$logLik, n*log(n/areas$totalArea) - n)
})

test_that("Scaling Z does not change the fit with baseline", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$withZeros
  fit <- fit.given.covar.baseline(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                  returnParams = TRUE)
  # With b around 1, (1e300)^b overflows and (1e-300)^b underflows
  for(scale in c(1e300, 1e-300)){
    fitScaled <- fit.given.covar.baseline(Z = scale*Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                          returnParams = TRUE)
    expect_equal(fitScaled$b, fit$b, tolerance = 1e-4)
    expect_equal(fitScaled$w, fit$w, tolerance = 1e-4)
    expect_equal(fitScaled$logLik, fit$logLik, tolerance = 1e-6)
  }
})

test_that("Invalid input to the model with baseline gives an error", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  Z <- make.baseline.test.covariates(setup = setup, tidspunkt = tidspunkt)$withZeros
  call.with <- function(Z, ...){
    fit.given.covar.baseline(Z = Z, structures = points$scars, l_x = setup$l_x, l_y = setup$l_y, returnParams = FALSE,
                             ...)
  }

  expect_error(call.with(Z = Z - 0.1), "non-negative")
  expect_error(call.with(Z = 0*Z), "zero everywhere")
  expect_error(call.with(Z = Z, bGrid = c(-1, 0, 1)), "at least 0")
})
