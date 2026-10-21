# Log likelihood from mppm with one dummy point per pixel centre, which makes the spatstat likelihood exact
exact.spatstat.logLik <- function(setup, tidspunkt, h, vindauga, skall = NULL, localStat = "variance"){
  field_list <- lapply(X = tidspunkt + h, FUN = surf.transform.specific.t, field = setup$field, l_x = setup$l_x,
                       l_y = setup$l_y, vindauga = vindauga, skall = skall, localStat = localStat)
  to.quad <- function(points){
    pp <- make.points.list.spatstat(points_list = points, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
    mapply(function(X, Z){spatstat.geom::pixelquad(X, W = spatstat.geom::as.owin(Z))}, pp, field_list, SIMPLIFY = FALSE)
  }
  df <- spatstat.geom::hyperframe(scars = to.quad(setup$scars), dimples = to.quad(setup$dimples), field = field_list)
  return(c(scars = as.numeric(stats::logLik(spatstat.model::mppm(scars ~ log(field), data = df))),
           dimples = as.numeric(stats::logLik(spatstat.model::mppm(dimples ~ log(field), data = df)))))
}

test_that("The local model matches the exact spatstat likelihood", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  vindauga <- t(finn.vindauga.sirkel(radius = 2))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for(localStat in c("variance", "meanSquare", "meanAbs")){
    res <- get.logLik.specific.r.h(r = 2, h = 1, field = setup$field, scars = points$scars, dimples = points$dimples,
                                   l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt, localStat = localStat)
    expected <- exact.spatstat.logLik(setup = setup, tidspunkt = tidspunkt, h = 1, vindauga = vindauga, skall = skall,
                                      localStat = localStat)
    expect_equal(res, expected, tolerance = 1e-6)
  }
})

test_that("The global model matches the exact spatstat likelihood", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)

  for(localStat in c("variance", "meanSquare", "meanAbs")){
    res <- get.logLik.specific.r.h(r = "infty", h = -1, field = setup$field, scars = points$scars,
                                   dimples = points$dimples, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt,
                                   localStat = localStat)
    expected <- exact.spatstat.logLik(setup = setup, tidspunkt = tidspunkt, h = -1, vindauga = "infty",
                                      localStat = localStat)
    expect_equal(res, expected, tolerance = 1e-6)
  }
})

test_that("fit.global.covar gives the same fit as fit.given.covar with a constant covariate", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  # The snapshots 1:3, where the estimate of b for the dimples is about 42, so the search interval is also widened
  Zglobal <- vapply(tidspunkt - 1, function(t){mean(setup$field[t, , ]^2)}, numeric(1))
  Z <- array(data = rep(Zglobal, times = dim(setup$field)[2]*dim(setup$field)[3]),
             dim = c(length(tidspunkt), dim(setup$field)[2], dim(setup$field)[3]))

  for(structures in points){
    fitGlobal <- fit.global.covar(Zglobal = Zglobal, structures = structures, l_x = setup$l_x, l_y = setup$l_y,
                                  returnParams = TRUE)
    fitFull <- fit.given.covar(Z = Z, structures = structures, l_x = setup$l_x, l_y = setup$l_y, returnParams = TRUE)
    expect_equal(fitGlobal$b, fitFull$b, tolerance = 1e-4)
    expect_equal(fitGlobal$a, fitFull$a, tolerance = 1e-3)
    expect_equal(fitGlobal$logLik, fitFull$logLik, tolerance = 1e-6)
  }

  expect_error(fit.global.covar(Zglobal = Zglobal[1:2], structures = points$scars, l_x = setup$l_x, l_y = setup$l_y,
                                returnParams = FALSE), "one element per snapshot")
})

test_that("make.global.result.vectors has the log likelihood for each lag", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  lags <- c(-1, 0, 2)
  res <- make.global.result.vectors(lags = lags, field = setup$field, scars = setup$scars, dimples = setup$dimples,
                                    l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
  expect_named(res, c("scars", "dimples"))
  expect_named(res$scars, c("-1", "0", "2"))

  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  for(h_ind in seq_along(lags)){
    expected <- get.logLik.global.h(h = lags[h_ind], field = setup$field, scars = points$scars, dimples = points$dimples,
                                    l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
    expect_equal(c(scars = res$scars[[h_ind]], dimples = res$dimples[[h_ind]]), expected)
  }

  single <- make.global.result.vectors(lags = 0, field = setup$field, scars = setup$scars, dimples = setup$dimples,
                                       l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
  expect_length(single$dimples, 1)
})

test_that("fit.given.covar works when Z^b overflows or underflows", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  # Global covariate for the snapshots 1:3, where the estimate of b for the dimples is about 42
  field <- setup$field[tidspunkt - 1, , , drop = FALSE]
  Z <- array(data = NA_real_, dim = dim(field))
  for(i in 1:length(tidspunkt)){
    Z[i, , ] <- mean(field[i, , ]^2)
  }
  fit <- fit.given.covar(Z = Z, structures = points$dimples, l_x = setup$l_x, l_y = setup$l_y, returnParams = TRUE)
  expect_gt(fit$b, 20)

  # Scaling Z only changes a, and with b around 42, (1e8)^b overflows and (1e-8)^b underflows
  for(scale in c(1e8, 1e-8)){
    fitScaled <- fit.given.covar(Z = scale*Z, structures = points$dimples, l_x = setup$l_x, l_y = setup$l_y,
                                 returnParams = TRUE)
    expect_equal(fitScaled$b, fit$b, tolerance = 1e-4)
    expect_equal(fitScaled$logLik, fit$logLik, tolerance = 1e-6)
  }

  Z[1, 1, 1] <- 0
  expect_error(fit.given.covar(Z = Z, structures = points$dimples, l_x = setup$l_x, l_y = setup$l_y, returnParams = FALSE),
               "must be positive")
})

test_that("Radius and lag encoded in one number give the same result", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  encoded <- make.input.matrix(radiar = 3, lags = 2)[1, 1]

  expect_equal(get.logLik.specific.r.h(r = encoded, h = NULL, field = setup$field, scars = points$scars,
                                       dimples = points$dimples, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt),
               get.logLik.specific.r.h(r = 3, h = 2, field = setup$field, scars = points$scars,
                                       dimples = points$dimples, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt))
})

test_that("A single snapshot keeps its time dimension", {
  setup <- make.likelihood.test.setup()
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = 2)

  res <- get.logLik.specific.r.h(r = 2, h = 0, field = setup$field, scars = points$scars, dimples = points$dimples,
                                 l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 2)
  expect_named(res, c("scars", "dimples"))
  expect_true(all(is.finite(res)))
})

test_that("make.result.matrices has the log likelihood for each lag and radius", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  radiar <- c(1, 2, 3)
  lags <- c(-1, 0, 2)
  res <- make.result.matrices(radiar = radiar, lags = lags, field = setup$field, scars = setup$scars,
                              dimples = setup$dimples, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
  expect_named(res, c("scars", "dimples"))
  expect_equal(dim(res$scars), c(length(lags), length(radiar)))
  expect_equal(dim(res$dimples), c(length(lags), length(radiar)))

  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  for(h_ind in seq_along(lags)){
    for(r_ind in seq_along(radiar)){
      expected <- get.logLik.specific.r.h(r = radiar[r_ind], h = lags[h_ind], field = setup$field, scars = points$scars,
                                          dimples = points$dimples, l_x = setup$l_x, l_y = setup$l_y,
                                          tidspunkt = tidspunkt)
      expect_equal(c(scars = res$scars[h_ind, r_ind], dimples = res$dimples[h_ind, r_ind]), expected)
    }
  }
})

test_that("make.result.matrices returns matrices for a single lag or radius", {
  setup <- make.likelihood.test.setup()
  call.with <- function(radiar, lags){
    make.result.matrices(radiar = radiar, lags = lags, field = setup$field, scars = setup$scars, dimples = setup$dimples,
                         l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 2:4)
  }

  expect_equal(dim(call.with(radiar = c(1, 2), lags = 0)$dimples), c(1, 2))
  res <- call.with(radiar = 2, lags = c(0, 1))
  expect_equal(dim(res$scars), c(2, 1))
  expect_equal(dimnames(res$scars), list(lag = c("0", "1"), radius = "2"))
  expect_error(call.with(radiar = 1.5, lags = 0), "whole numbers")
})

test_that("Invalid input gives an error", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)
  call.with <- function(...){
    get.logLik.specific.r.h(field = setup$field, scars = points$scars, dimples = points$dimples,
                            l_x = setup$l_x, l_y = setup$l_y, ...)
  }

  expect_error(call.with(r = 2, h = 0, tidspunkt = tidspunkt, localStat = "mean"))
  expect_error(call.with(r = 0, h = 0, tidspunkt = tidspunkt), "at least one pixel")
  expect_error(call.with(r = 2, h = 0, tidspunkt = 2:3), "one element per timestep")
})

test_that("Empty timesteps can be stored as 0 or as a data frame with zero rows", {
  pointsList <- list(0, structures.to.data.frame(list()), data.frame(x = 0.5, y = 0.25))
  res <- make.points.list(pointsList = pointsList, l_x = 2, l_y = 4, tidspunkt = 1:3)
  expect_equal(sapply(res, nrow), c(0, 0, 1))
  expect_equal(res[[3]], data.frame(x = 1, y = 1))

  resSpatstat <- make.points.list.spatstat(points_list = pointsList, l_x = 2, l_y = 4, tidspunkt = 1:3)
  expect_equal(sapply(resSpatstat, spatstat.geom::npoints), c(0, 0, 1))
})
