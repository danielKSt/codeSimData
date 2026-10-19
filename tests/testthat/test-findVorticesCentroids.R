test_that("Simple Vortex is explored correctly", {
  simplePixels <- data.frame(x = c(2,3), y = c(2,2))
  expect_equal(explore.vortex(vortexPixels = simplePixels, boxDimensions = c(4, 3)),
              c((2.5-0.5)/4, (2-0.5)/3, 2, 1, 0))
})

test_that("Vortex crossing border is explored correctly", {
  crossingPixels <- data.frame(x = c(1, 1, 5), y = c(2, 3, 2))
  expect_equal(explore.vortex(vortexPixels = crossingPixels, boxDimensions = c(5, 5)),
               c((2/3-0.5)/5, (2+1/3-0.5)/5, 3, 1, 1))
})

test_that("Vortex crossing border with centre in the last pixel wraps around", {
  crossingPixels <- data.frame(x = c(1, 5, 5), y = c(2, 2, 3))
  expect_equal(explore.vortex(vortexPixels = crossingPixels, boxDimensions = c(5, 5)),
               c(1-1/30, (2+1/3-0.5)/5, 3, 1, 1))
})

test_that("Scar centres are placed at pixel centres", {
  # The skeleton of a line one pixel wide is the line itself, so we skip the call to python
  local_mocked_bindings(skeletonize_object_py = function(df, np, skimage) df)

  scarPixels <- data.frame(x = c(3, 4, 5), y = c(2, 2, 2))
  expect_equal(explore.scar(pixels = scarPixels, dims = c(8, 6), np = "unused", skimage = "unused"),
               c(3.5/8, 1.5/6, 3, 2, 0))

  crossingPixels <- data.frame(x = c(7, 8, 1), y = c(3, 3, 3))
  expect_equal(explore.scar(pixels = crossingPixels, dims = c(8, 6), np = "unused", skimage = "unused"),
               c(7.5/8, 2.5/6, 3, 2, 0))

  crossingPixelsOddDim <- data.frame(x = c(6, 7, 1), y = c(3, 3, 3))
  expect_equal(explore.scar(pixels = crossingPixelsOddDim, dims = c(7, 5), np = "unused", skimage = "unused"),
               c(6.5/7, 2.5/5, 3, 2, 0))
})
