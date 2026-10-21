# The field is not square, so mixing up the x and y directions makes the tests fail
make.likelihood.test.setup <- function(){
  set.seed(11)
  nx <- 12
  ny <- 10
  noTimesteps <- 6
  field <- array(data = stats::rnorm(noTimesteps*nx*ny), dim = c(noTimesteps, nx, ny))
  random.points <- function(n){data.frame(x = stats::runif(n), y = stats::runif(n))}
  scars <- lapply(1:noTimesteps, function(t){random.points(10)})
  dimples <- lapply(1:noTimesteps, function(t){random.points(8)})
  # An empty timestep, stored the same way as by find.vortices.over.time
  dimples[[3]] <- 0
  return(list(field = field, scars = scars, dimples = dimples, l_x = nx/2.5, l_y = ny/2.5))
}

prepare.likelihood.test.points <- function(setup, tidspunkt){
  list(scars = make.points.list(pointsList = setup$scars, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt),
       dimples = make.points.list(pointsList = setup$dimples, l_x = setup$l_x, l_y = setup$l_y, tidspunkt = tidspunkt))
}
