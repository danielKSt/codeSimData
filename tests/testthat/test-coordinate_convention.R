# A structure covering a single pixel should, after being turned into a point, look up the covariate
# value of that same pixel. We use a non-square field and a side length different from one, so that
# swapped axes or a wrong pixel size show up.

make.test.setup <- function(){
  nx <- 8
  ny <- 6
  noTimesteps <- 3
  field <- array(data = 0, dim = c(noTimesteps, nx, ny))
  for(t in 1:noTimesteps){
    field[t, , ] <- outer(1:nx, 1:ny, function(i, j){100*t + 10*i + j})
  }
  structurePixels <- list(rbind(c(1, 1), c(5, 3)),
                          rbind(c(8, 6), c(3, 5)),
                          rbind(c(1, 6), c(8, 3)))
  mask <- array(data = 0L, dim = c(noTimesteps, nx, ny))
  for(t in 1:noTimesteps){
    mask[t, , ][structurePixels[[t]]] <- 1L
  }
  return(list(field = field, mask = mask, structurePixels = structurePixels, l_x = nx/2.5, l_y = ny/2.5))
}

test_that("Vortex points look up the covariate value of their own pixel", {
  setup <- make.test.setup()
  vortices <- find.vortices.over.time(vortexData = setup$mask)
  points <- make.points.list.spatstat(points_list = vortices, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 1:3)

  for(t in 1:3){
    Z <- surf.transform.specific.t(t = t, field = setup$field, l_x = setup$l_x, l_y = setup$l_y)
    expect_equal(sort(Z[points[[t]]]), sort(setup$field[t, , ][setup$structurePixels[[t]]]))
  }
})

test_that("Scar points look up the covariate value of their own pixel", {
  # The skeleton of a single pixel is the pixel itself, so we skip the call to python
  local_mocked_bindings(skeletonize_object_py = function(df, np, skimage) df)
  setup <- make.test.setup()
  scars <- lapply(1:3, function(t){
    find.scars.specific.time(scarData = setup$mask[t, , ], np = "unused", skimage = "unused")
  })
  points <- make.points.list.spatstat(points_list = scars, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 1:3)

  for(t in 1:3){
    Z <- surf.transform.specific.t(t = t, field = setup$field, l_x = setup$l_x, l_y = setup$l_y)
    expect_equal(sort(Z[points[[t]]]), sort(setup$field[t, , ][setup$structurePixels[[t]]]))
  }
})

test_that("Points look up the local variance of their own pixel", {
  setup <- make.test.setup()
  vortices <- find.vortices.over.time(vortexData = setup$mask)
  points <- make.points.list.spatstat(points_list = vortices, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 1:3)
  vindauga <- t(finn.vindauga.sirkel(radius = 2))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for(t in 1:3){
    Z <- surf.transform.specific.t(t = t, field = setup$field, l_x = setup$l_x, l_y = setup$l_y,
                                   vindauga = vindauga, skall = skall)
    locVar <- locvar.transformation.periodic(sField = setup$field[t, , ], vindauga = vindauga, skall = skall)
    expect_equal(sort(Z[points[[t]]]), sort(locVar[setup$structurePixels[[t]]]))
  }
})

test_that("Points in the stacked window look up the covariate value of their own pixel", {
  withr::local_options(codeSimData.silenceDeprecation = TRUE)
  setup <- make.test.setup()
  vortices <- find.vortices.over.time(vortexData = setup$mask)
  points <- prepare.points.dns(indices = 1:3, points = vortices)
  Z <- combine.covariates.dns(hDiv = setup$field, indices = 1:3, r = 2, h = 0)
  vindauga <- t(finn.vindauga.sirkel(radius = 2))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  expected <- unlist(lapply(1:3, function(t){
    locVar <- locvar.transformation.periodic(sField = setup$field[t, , ], vindauga = vindauga, skall = skall)
    return(locVar[setup$structurePixels[[t]]])
  }))
  expect_equal(sort(Z[points]), sort(expected))
})

test_that("get_Z_at_coord looks up the covariate value of the point's own pixel", {
  setup <- make.test.setup()
  vortices <- find.vortices.over.time(vortexData = setup$mask)
  vortices[[3]] <- 0
  points <- make.points.list(pointsList = vortices, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = 1:3)
  # An empty last timestep must not be dropped from the list
  expect_length(points, 3)
  expect_equal(nrow(points[[3]]), 0)

  expected <- list()
  for(t in 1:2){
    Zvals <- sapply(seq_len(nrow(points[[t]])), function(i){
      get_Z_at_coord(Z = setup$field, coord = c(points[[t]]$x[i], points[[t]]$y[i]), t = t,
                     l_x = setup$l_x, l_y = setup$l_y)
    })
    expected[[t]] <- setup$field[t, , ][setup$structurePixels[[t]]]
    expect_equal(sort(Zvals), sort(expected[[t]]))
  }

  inData <- get_sumLogZ_and_n(Z = setup$field, structures = points, l_x = setup$l_x, l_y = setup$l_y)
  expect_equal(inData$n, 4)
  expect_equal(inData$sumLogZ, sum(log(unlist(expected))))

  # The far edges of the window belong to the last pixel, and points outside the window are an error
  expect_equal(get_Z_at_coord(Z = setup$field, coord = c(setup$l_x, setup$l_y), t = 1, l_x = setup$l_x, l_y = setup$l_y),
               setup$field[1, 8, 6])
  expect_error(get_Z_at_coord(Z = setup$field, coord = c(-0.1, 1), t = 1, l_x = setup$l_x, l_y = setup$l_y),
               "outside the observation window")
})

test_that("Old vortex and scar coordinates are converted to pixel centres", {
  dims <- c(8, 6)
  oldVortices <- list(data.frame(x = c(1, 8)/(8+1), y = c(6, 2.5)/(6+1), noPixels = c(1, 2)),
                      0,
                      data.frame(x = 1/3/(8+1), y = 3/(6+1), noPixels = 3))
  newVortices <- convert.points.to.pixel.centres(points_list = oldVortices, dims = dims, type = "vortices")
  expect_equal(newVortices[[1]]$x, c(0.5, 7.5)/8)
  expect_equal(newVortices[[1]]$y, c(5.5, 2)/6)
  expect_equal(newVortices[[1]]$noPixels, c(1, 2))
  expect_equal(newVortices[[2]], 0)
  expect_equal(newVortices[[3]]$x, (1/3-0.5)/8 + 1)

  oldScars <- list(list(data.frame(x = (c(1, 8)/8) %% 1, y = (c(6, 3)/6) %% 1)))
  newScars <- convert.points.to.pixel.centres(points_list = oldScars, dims = dims, type = "scars")
  expect_equal(newScars[[1]][[1]]$x, c(0.5, 7.5)/8)
  expect_equal(newScars[[1]][[1]]$y, c(5.5, 2.5)/6)

  expect_error(convert.points.to.pixel.centres(points_list = newVortices, dims = dims, type = "vortices"),
               "already use the pixel-centre convention")

  scaledTwice <- list(data.frame(x = ((c(1, 8)/8) %% 1)/8, y = c(0.5, 0.5)))
  expect_warning(convert.points.to.pixel.centres(points_list = scaledTwice, dims = dims, type = "scars"),
                 "not on the old pixel grid")
})

test_that("Scars crossing the border that were divided twice are repaired before conversion", {
  dims <- c(8, 6)
  # Old explore.scar returned unit coordinates in the directions a scar crossed the border, and pixel indices
  # otherwise, and everything was then divided by the number of pixels.
  # Scar in pixel (7, 6) crossing both borders: (7/8, 0) divided again.
  # Scar in pixel (8, 2) crossing no border: (8, 2) divided once.
  # Scar in pixel (1, 4) crossing the x border only: (1/8, 4) divided once more in x and once in y.
  dividedTwice <- list(data.frame(x = c(7/64, 1, 1/64), y = c(0, 2/6, 4/6)), 0)

  expect_no_warning(repaired <- convert.points.to.pixel.centres(points_list = dividedTwice, dims = dims,
                                                                type = "scars", repairSeam = TRUE))
  expect_equal(repaired[[1]]$x, c(6.5, 7.5, 0.5)/8)
  expect_equal(repaired[[1]]$y, c(5.5, 1.5, 3.5)/6)
  expect_equal(repaired[[2]], 0)

  expect_warning(convert.points.to.pixel.centres(points_list = dividedTwice, dims = dims, type = "scars"),
                 "repairSeam = TRUE")
  expect_error(convert.points.to.pixel.centres(points_list = dividedTwice, dims = dims, type = "vortices",
                                               repairSeam = TRUE),
               "only be used")
})
