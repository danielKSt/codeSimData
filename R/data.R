#' Numerically Simulated Surface Divergence data
#' Generated using the command:
#'    library("rhdf5")
#'    surfaceDivergence <- h5read(file = '/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/surface_vortices.h5', name = 'sdiv')
#'
#' @format ## surfaceDivergence
#' A large 3-dimensional array with 12500 timesteps of surface divergence in a 256-by-256 grid for every timestep
"surfaceDivergence"

#' Numerically Simulated Surface Vortex data
#' Generated using the commands:
#'    library("rhdf5")
#'    surfaceVortices <- h5read(file = '/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/surface_vortices.h5', name = 'surf')
#'
#' @format ## surfaceDivergence
#' A large 3-dimensional array with 12500 timesteps in a 256-by-256 grid for every timestep every pixel is either a one or a zero, with a one if there is a vortex, and a zero if there is not.
"surfaceVortices"


# surfaceDivergence <- h5read(file = '/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/surface_vortices.h5', name = 'sdiv')
# surfaceVortices <- h5read(file = '/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/SimData/data/surface_vortices.h5', name = 'surf')
# surfaceDivergence <- surfaceDivergence[1:5, , ]
# surfaceVortices <- surfaceVortices[1:5, , ]
# usethis::use_data(surfaceVortices, surfaceDivergence, overwrite = TRUE)
