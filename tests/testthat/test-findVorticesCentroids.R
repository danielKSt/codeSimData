test_that("Simple Vortex is explored correctly", {
  simplePixels <- data.frame(x = c(2,3), y = c(2,2))
  expect_equal(explore.vortex(vortexPixels = simplePixels, boxDimensions = c(4, 3)),
              c(2.5/5, 2/4, 2, 1, 0))
})

test_that("Vortex crossing border is explored correctly", {
  crossingPixels <- data.frame(x = c(1, 1, 5), y = c(2, 3, 2))
  expect_equal(explore.vortex(vortexPixels = crossingPixels, boxDimensions = c(5, 5)),
               c(2/(3*(5+1)), (2+1/3)/(5+1), 3, 1, 1))
})
