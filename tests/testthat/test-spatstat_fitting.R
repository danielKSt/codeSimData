# With one dummy point per pixel centre, mppm gives the exact likelihood, which is what result_matrices.R calculates

test_that("get.fit.specific.r.h with pixel quadrature matches the exact likelihood", {
  setup <- make.likelihood.test.setup()
  tidspunkt <- 2:4
  points <- prepare.likelihood.test.points(setup = setup, tidspunkt = tidspunkt)

  for(r in list(2, "infty")){
    expected <- get.logLik.specific.r.h(r = r, h = 1, field = setup$field, scars = points$scars, dimples = points$dimples,
                                        l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt)
    res <- get.fit.specific.r.h(scars_list = setup$scars, dimples_list = setup$dimples, l_x = setup$l_x, l_y = setup$l_y,
                                field = setup$field, tidspunkt = tidspunkt, r = r, h = 1)
    expect_equal(res, expected, tolerance = 1e-6)
  }
})

test_that("make.logLik.matrix matches make.result.matrices", {
  setup <- make.likelihood.test.setup()
  args <- list(field = setup$field, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 2:4, radiar = c(1, 2), lags = c(0, 1))

  res <- do.call(make.logLik.matrix, c(args, list(scars_list = setup$scars, dimples_list = setup$dimples,
                                                  radiiForPrint = -1)))
  expected <- do.call(make.result.matrices, c(args, list(scars = setup$scars, dimples = setup$dimples)))
  expect_equal(res, expected, tolerance = 1e-6)

  # A single lag still gives a lag x radius matrix
  args$lags <- 1
  res <- do.call(make.logLik.matrix, c(args, list(scars_list = setup$scars, dimples_list = setup$dimples,
                                                  radiiForPrint = -1, pixelQuadrature = FALSE)))
  expect_equal(dim(res$scars), c(1, 2))
  expect_equal(dimnames(res$scars), list(lag = "1", radius = c("1", "2")))
})

test_that("get.fit.specific.r.h names its results", {
  setup <- make.likelihood.test.setup()
  call.with <- function(...){
    get.fit.specific.r.h(scars_list = setup$scars, dimples_list = setup$dimples, l_x = setup$l_x, l_y = setup$l_y,
                         field = setup$field, tidspunkt = 2:4, r = 2, h = 0, pixelQuadrature = FALSE, ...)
  }

  expect_named(call.with(), c("scars", "dimples"))
  expect_named(call.with(fitLogless = TRUE), c("scars", "dimples", "scars.noLog", "dimples.noLog"))
  expect_named(call.with(resType = "aic"), c("scars", "dimples"))
})

test_that("Invalid input to the spatstat fitting is rejected", {
  setup <- make.likelihood.test.setup()

  # The raw field is signed, so log(field) is undefined
  expect_error(get.fit.specific.r.h(scars_list = setup$scars, dimples_list = setup$dimples, l_x = setup$l_x,
                                    l_y = setup$l_y, field = setup$field, tidspunkt = 2:4, r = 0, h = 0),
               "must be positive")
  expect_error(make.logLik.matrix(scars_list = setup$scars, dimples_list = setup$dimples, l_x = setup$l_x,
                                  l_y = setup$l_y, field = setup$field, tidspunkt = 2:4, radiar = 1.5, lags = 0,
                                  radiiForPrint = -1),
               "whole numbers")
  expect_error(make.input.matrix(radiar = 1, lags = 0.5), "whole numbers")
})
