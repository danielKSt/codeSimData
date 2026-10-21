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

test_that("Connected structures are found with and without periodic boundaries", {
  mask <- matrix(0L, nrow = 6, ncol = 5)
  # Neighbours across the periodic border in x, and diagonal neighbours
  mask[cbind(c(1, 6, 3, 4), c(2, 2, 3, 4))] <- 1L
  sorted.pixels <- function(component){as.matrix(component[order(component$x, component$y), ])}

  periodic <- find.connected.components(mask = mask, periodic = TRUE)
  expect_length(periodic, 2)
  expect_equal(sorted.pixels(periodic[[1]]), cbind(x = c(1, 6), y = c(2, 2)), ignore_attr = TRUE)
  expect_equal(sorted.pixels(periodic[[2]]), cbind(x = c(3, 4), y = c(3, 4)), ignore_attr = TRUE)

  # Without periodic boundaries the pixels on opposite edges are separate, and the structures are ordered by their
  # first pixel with x in the outer loop
  physical <- find.connected.components(mask = mask, periodic = FALSE)
  expect_length(physical, 3)
  expect_equal(sapply(physical, function(component){min(component$x)}), c(1, 3, 6))
})

test_that("A structure in the same x index as an earlier structure is found", {
  # Earlier versions continued scanning with the x index of the last pixel of the previous structure, and missed
  # the pixel (2, 4)
  mask <- matrix(0L, nrow = 4, ncol = 6)
  mask[cbind(c(2, 3, 2), c(1, 1, 4))] <- 1L
  vortices <- find.vortices.specific.time(vortexData = mask, periodic = FALSE)
  expect_equal(vortices$noPixels, c(2, 1))
})

test_that("Timesteps without structures give data frames with zero rows", {
  vortices <- find.vortices.specific.time(vortexData = matrix(0L, nrow = 4, ncol = 4))
  expect_equal(nrow(vortices), 0)
  expect_named(vortices, c("x", "y", "noPixels", "xStretch", "yStretch"))

  expect_true(is.empty.snapshot(vortices))
  expect_true(is.empty.snapshot(0))
  expect_false(is.empty.snapshot(data.frame(x = 0.5, y = 0.5)))
  expect_false(is.empty.snapshot(list(vortices)))
})

test_that("Scars are found from the mask", {
  # The skeleton of a line one pixel wide is the line itself, so we skip the calls to python
  local_mocked_bindings(skeletonize_object_py = function(df, np, skimage) df,
                        python.modules = function(np = NULL, skimage = NULL) list(np = "unused", skimage = "unused"))
  mask <- matrix(0L, nrow = 8, ncol = 6)
  mask[3:5, 2] <- 1L

  expect_equal(find.scars.specific.time(scarData = mask), structures.to.data.frame(list(c(3.5/8, 1.5/6, 3, 2, 0))))
  expect_length(find.scars.centrelines.specific.time(scarData = mask), 1)
  expect_equal(nrow(find.scars.specific.time(scarData = 0L*mask)), 0)
  expect_length(find.scars.centrelines.specific.time(scarData = 0L*mask), 0)
})
