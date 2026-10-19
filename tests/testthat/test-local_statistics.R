# The sliding window calculations in C++ should give the same as calculating each window directly. We use a
# non-square field with non-zero mean, so that swapped axes, or mixing up the mean square and the variance, show up.
# The field also has negative values, so that the mean of absolute values differs from the mean.

make.test.field <- function(){
  set.seed(1)
  return(matrix(data = stats::rnorm(13*11, mean = 0.5), nrow = 13, ncol = 11))
}

brute.force.local.statistic <- function(sField, radius, statistic, periodic){
  nx <- dim(sField)[1]
  ny <- dim(sField)[2]
  vindauga <- finn.vindauga.sirkel(radius = radius)
  res <- matrix(data = NA, nrow = nx, ncol = ny)
  for(i in 1:nx){
    for(j in 1:ny){
      xi <- i + vindauga[, 1]
      yj <- j + vindauga[, 2]
      if(periodic){
        xi <- (xi - 1) %% nx + 1
        yj <- (yj - 1) %% ny + 1
      } else {
        inside <- (xi >= 1) & (xi <= nx) & (yj >= 1) & (yj <= ny)
        xi <- xi[inside]
        yj <- yj[inside]
      }
      res[i, j] <- statistic(sField[cbind(xi, yj)])
    }
  }
  return(res)
}

sliding.local.statistic <- function(calculateLocStat, sField, radius){
  vindauga <- t(finn.vindauga.sirkel(radius = radius))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))
  res <- calculateLocStat(inWin = vindauga, inShell = skall, insDiv = sField,
                          inDims = c(dim(sField)[1], dim(sField)[2], dim(vindauga)[2], dim(skall)[2]))
  return(matrix(data = res, nrow = dim(sField)[1]))
}

statistics <- list(
  variance = list(brute = stats::var, periodic = calculateLocVar_periodic, physical = calculateLocVar_physical),
  meanSquare = list(brute = function(x){mean(x^2)}, periodic = calculateLocMeanSq_periodic,
                    physical = calculateLocMeanSq_physical),
  meanAbs = list(brute = function(x){mean(abs(x))}, periodic = calculateLocMeanAbs_periodic,
                 physical = calculateLocMeanAbs_physical),
  mean = list(brute = mean, periodic = calculateLocMean_periodic, physical = calculateLocMean_physical)
)

test_that("Local statistics with periodic boundary match direct calculation", {
  sField <- make.test.field()
  for(statName in names(statistics)){
    for(radius in 1:5){
      expect_equal(sliding.local.statistic(statistics[[statName]]$periodic, sField = sField, radius = radius),
                   brute.force.local.statistic(sField = sField, radius = radius,
                                               statistic = statistics[[statName]]$brute, periodic = TRUE),
                   label = paste(statName, "with radius", radius))
    }
  }
})

test_that("Local statistics with physical boundary match direct calculation", {
  sField <- make.test.field()
  for(statName in names(statistics)){
    for(radius in 1:5){
      expect_equal(sliding.local.statistic(statistics[[statName]]$physical, sField = sField, radius = radius),
                   brute.force.local.statistic(sField = sField, radius = radius,
                                               statistic = statistics[[statName]]$brute, periodic = FALSE),
                   label = paste(statName, "with radius", radius))
    }
  }
})

test_that("surf.transform.specific.t uses the chosen local statistic", {
  sField <- make.test.field()
  field <- array(data = 0, dim = c(2, dim(sField)))
  field[2, , ] <- sField
  fieldEnsemble <- array(data = 0, dim = c(2, 2, dim(sField)))
  fieldEnsemble[1, 2, , ] <- sField
  vindauga <- t(finn.vindauga.sirkel(radius = 3))
  skall <- t(finn.skall.av.vindauga(vindauga = t(vindauga)))

  for(statName in names(statistics)){
    Z <- surf.transform.specific.t(t = 2, field = field, l_x = 2, l_y = 1.5, vindauga = vindauga, skall = skall,
                                   localStat = statName)
    expected <- brute.force.local.statistic(sField = sField, radius = 3,
                                            statistic = statistics[[statName]]$brute, periodic = TRUE)
    expect_equal(Z$v, t(expected), ignore_attr = TRUE, label = paste(statName, "periodic"))

    Z <- surf.transform.specific.t(t = c(1, 2), field = fieldEnsemble, l_x = 2, l_y = 1.5, vindauga = vindauga,
                                   skall = skall, localStat = statName)
    expected <- brute.force.local.statistic(sField = sField, radius = 3,
                                            statistic = statistics[[statName]]$brute, periodic = FALSE)
    expect_equal(Z$v, t(expected), ignore_attr = TRUE, label = paste(statName, "physical"))
  }

  Z <- surf.transform.specific.t(t = 2, field = field, l_x = 2, l_y = 1.5, vindauga = "infty", localStat = "mean")
  expect_equal(unique(as.vector(Z$v)), mean(sField))
  Z <- surf.transform.specific.t(t = 2, field = field, l_x = 2, l_y = 1.5, vindauga = "infty", localStat = "meanSquare")
  expect_equal(unique(as.vector(Z$v)), mean(sField^2))
  Z <- surf.transform.specific.t(t = 2, field = field, l_x = 2, l_y = 1.5, vindauga = "infty", localStat = "meanAbs")
  expect_equal(unique(as.vector(Z$v)), mean(abs(sField)))

  expect_error(surf.transform.specific.t(t = 2, field = field, l_x = 2, l_y = 1.5, vindauga = vindauga, skall = skall,
                                         localStat = "median"))
})
