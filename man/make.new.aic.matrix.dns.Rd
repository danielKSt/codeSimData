% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aic_matrices.R
\name{make.new.aic.matrix.dns}
\alias{make.new.aic.matrix.dns}
\title{Function for making a new AIC matrix given DNS data}
\usage{
make.new.aic.matrix.dns(
  radiar,
  lags,
  vortices.pp,
  scars.pp,
  sDiv,
  indices,
  radiiModulo = 1000L,
  radiiForPrint = c(1)
)
}
\arguments{
\item{radiar}{Vector with all radii for local variance}

\item{lags}{Vector with all lag values to use}

\item{vortices.pp}{votrices as spatstat point process}

\item{scars.pp}{scars as spatstat point process}

\item{sDiv}{horisontal divergence}

\item{indices}{which indices to use}

\item{radiiModulo}{parameter for encoding lag and radius into single matrix}

\item{radiiForPrint}{Boolean variable to indicate whether to print progress or not}
}
\description{
Function for making a new AIC matrix given DNS data
}
