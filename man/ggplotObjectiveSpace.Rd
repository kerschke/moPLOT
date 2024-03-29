% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplotObjectiveSpace.R
\name{ggplotObjectiveSpace}
\alias{ggplotObjectiveSpace}
\title{Visualize Multi-Objective Gradient Landscape in the Objective Space.}
\usage{
ggplotObjectiveSpace(
  df,
  var1 = "y1",
  var2 = "y2",
  log.scale = TRUE,
  impute.zero = TRUE,
  minimalistic.image = FALSE,
  color.palette,
  legend.position,
  ...
)
}
\arguments{
\item{df}{[\code{\link[=data.frame]{data.frame()}}]\cr
Data frame containing the positions of the individuals (per row) and their
corresponding cumulated path length (denoted: \code{"height"}) as returned by
\code{\link[=computeCumulatedPathLengths]{computeCumulatedPathLengths()}}.}

\item{var1}{[\verb{[character](1L)}]\cr
Name of the variable indicating the first objective (default: \code{"y1"}).}

\item{var2}{[\verb{[character](1L)}]\cr
Name of the variable indicating the second objective (default: \code{"y2"}).}

\item{log.scale}{[\verb{[logical](1L)}]\cr
Should the resulting heights be displayed on a log-scale? The default is \code{TRUE}.}

\item{impute.zero}{[\verb{[logical](1L)}]\cr
Should height values, which are exactly zero be imputed by a value half the magnitude
of the smallest non-zero height? Otherwise ggplot will automatically color the
corresponding tiles by a color representing \code{NA} values (usually grey).
Note that this parameter is only relevant in case of \code{log.scale = TRUE}.
The default is \code{TRUE}.}

\item{minimalistic.image}{[\verb{[logical](1L)}]\cr
Should all information surrounding the image (axes, legends, background, etc.) be discarded?
The default is \code{FALSE}.}

\item{color.palette}{[\code{\link[=character]{character()}}]\cr
Vector of colors used for visualizing the different heights of the landscape. By default,
this function tries to use the color palettes from \code{fields::tim.color} or
\code{viridis}. However, if neither of these packages is installed, it will use
\code{terrain.colors}.}

\item{legend.position}{[\verb{[character](1L)}]\cr
On which side of the plot, should the legend be located? If this information is not provided
and \code{minimalisitic.image = FALSE}, the legend will be placed on the right side.}

\item{...}{\link{any}\cr
Further arguments to be passed to the \code{geom_tile} function of \code{ggplot}.}
}
\value{
\code{\link{ggplot}}\cr
A \code{ggplot} object displaying the multi-objective gradient landscape.
}
\description{
Use the coloring of the cumulated path lengths (see
\code{\link[=computeCumulatedPathLengths]{computeCumulatedPathLengths()}}) and visualize the
corresponding points in the objective space.
}
