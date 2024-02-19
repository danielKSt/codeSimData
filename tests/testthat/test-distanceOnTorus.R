test_that("Distance is correct when d=1 0<x<y<0.5", {
  expect_equal(distanceOnTorus(1, 0.1, 0.2), 0.1)
})

test_that("Distance is correct when d=1 0<x<0.2<0.9<y", {
  expect_equal(distanceOnTorus(1, 0.1, 0.95), 0.15)
})

test_that("Distance is correct when d=1 0<y<0.2<0.9<x", {
  expect_equal(distanceOnTorus(1, 0.95, 0.1), 0.15)
})


test_that("Distance is correct when d=2, x=(0.2, 0.2), y=(0.4, 0.4)", {
  expect_equal(distanceOnTorus(d=2, x = c(0.2, 0.2), y = c(0.4,0.4)), sqrt(2)/5.0)
})

test_that("Distance is correct when d=2, x=(0.2, 0.2), y=(0.9, 0.8)", {
  expect_equal(distanceOnTorus(d=2, x = c(0.2, 0.2), y = c(0.9,0.8)), 0.50)
})

test_that("Distance is correct when d=2, x=(0.8, 0.9), y=(0.2, 0.2)", {
  expect_equal(distanceOnTorus(d=2, x = c(0.8, 0.9), y = c(0.2,0.2)), 0.50)
})
