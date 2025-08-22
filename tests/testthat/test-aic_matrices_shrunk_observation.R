test_that("Scaling into a shrunk window works", {
  lims <- c(0.25, 0.75)
  x_original <- c(0.1, 0.25, 0.375, 0.375, 0.5, 0.75, 0.9)
  y_original <- c(0.1, 0.25, 0.375, 0.75, 0.5, 0.75, 0.9)
  structures <- data.frame(x = x_original, y = y_original)

  inShrunk <- find.points.in.shrunk(points = structures, xlim = lims, ylim = lims)

  structures_shrunk <- scale.points.shrunk(points = structures[inShrunk, ], xlim = lims, ylim = lims)

  expect_equal(structures_shrunk$x, expected = c(0, 0.25, 0.25, 0.5, 1.0))
  expect_equal(structures_shrunk$y, expected = c(0, 0.25, 1.0, 0.5, 1.0))
})


test_that("Detecting structures outside of shrunk window works", {
  lims <- c(0.25, 0.75)
  x_original <- c(0.1, 0.25, 0.375, 0.375, 0.5, 0.75, 0.9)
  y_original <- c(0.1, 0.25, 0.375, 0.75, 0.5, 0.75, 0.9)
  structures <- data.frame(x = x_original, y = y_original)

  expect_equal(point.in.shrunk(points = structures, xlim = lims, ylim = lims), expected = TRUE)
  expect_equal(find.points.in.shrunk(points = structures, xlim = lims, ylim = lims), expected = c(2, 3, 4, 5, 6))

  expect_equal(point.in.shrunk(structures[c(1, 7), ], lims, lims), expected = FALSE)
})


test_that("Prepare structures for analysis DNS works", {
  lims <- c(0.25, 0.75)

  x_1 <- c(0.1, 0.25, 0.375, 0.375, 0.5, 0.75, 0.9)
  y_1 <- c(0.1, 0.25, 0.375, 0.75, 0.5, 0.75, 0.9)
  s_1 <- data.frame(x = x_1, y = y_1)

  x_2 <- c(0.1, 0.2, 0.8, 0.375, 0.1, 0.79, 0.9)
  y_2 <- c(0.1, 0.25, 0.375, 0.85, 0.5, 0.75, 0.9)
  s_2 <- data.frame(x = x_2, y = y_2)

  x_3 <- c(0.1, 0.65, 0.4, 0.6, 0.3, 0.76, 0.9)
  y_3 <- c(0.1, 0.25, 0.375, 0.75, 0.5, 0.75, 0.9)
  s_3 <- data.frame(x = x_3, y = y_3)

  structures <- list(s_1, s_2, s_3)

  res <- prepare.points.dns.shrunk(indices = c(1, 2, 3), points = structures, xlim = lims, ylim = lims)

  expect_equal(object = res$window$yrange, expected = c(0, 3))
  expect_equal(object = res$n, expected = 9)
  expect_equal(object = res$x, expected = c(0.0, 0.25, 0.25, 0.5, 1.0, 0.8, 0.3, 0.7, 0.1))
  expect_equal(object = res$y, expected = c(0.0, 0.25, 1.0, 0.5, 1.0, 2.0, 2.25, 3.0, 2.5))
})

# plot(x = structures$x, y = structures$y, xlim = c(0, 1), ylim = c(0, 1))
# lines(x = c(lims[1], lims[1], lims[2], lims[2], lims[1]), y = c(lims[1], lims[2], lims[2], lims[1], lims[1]), col = "red")
#
# plot(x = structures_shrunk$x, y = structures_shrunk$y, xlim = c(0, 1), ylim = c(0, 1))

test_that("Combine local variance for covariates", {
  sDiv_1 <- matrix(data = c(0, 1, 1, 0, 1,
                            1, 1, 0, 0, 1,
                            0, 2, 1, 1, 1,
                            1, 4, 0, 2, 1,
                            1, 0, 1, 0, 1), nrow = 5, ncol = 5)
  sDiv_2 <- matrix(data = c(0, 1, 1, 0, 1,
                            1, 1, 0, 0, 1,
                            0, 2, 1, 1, 1,
                            1, 4, 0, 2, 1,
                            1, 0, 1, 0, 1), nrow = 5, ncol = 5)

  sDiv <- array(data = c(sDiv_1, sDiv_2), dim = c(2, 5, 5))
  sDiv[1, , ] <- sDiv_1
  sDiv[2, , ] <- sDiv_2


  skim <- 1
  fullWindowSize <- 5

  lims <- c(skim + 1, fullWindowSize - skim)

  res <- combine.covariates.dns.shrunk(sDiv = sDiv, indices = c(1,2), r = 1, h = 0, xlim = lims, ylim = lims)

  expect_equal(object = res$v[, 1], expected = c(0.5, 2.3, 2.8, 0.5, 2.3, 2.8))
})





