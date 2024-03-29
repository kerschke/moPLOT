% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addGGOptima.R
\name{addGGOptima}
\alias{addGGOptima}
\title{Add (Known) Local Optima of Smoof Functions to Gradient Field Heatmap.}
\usage{
addGGOptima(g, fn1, fn2, fn3, symbol1 = 21L, symbol2 = 22L, symbol3 = 24L, ...)
}
\arguments{
\item{g}{[\code{\link[ggplot2:ggplot]{ggplot2::ggplot}}]\cr
\code{ggplot}-object as created by \code{\link{ggplotHeatmap}}.}

\item{fn1}{[\code{\link{function}}]\cr
The first objective used for computing the multi-objective gradient.}

\item{fn2}{[\code{\link{function}}]\cr
The second objective used for computing the multi-objective gradient.}

\item{fn3}{[\code{\link[=function]{function()}}]\cr
The third objective (if existing) used for computing the multi-objective gradient.}

\item{symbol1}{[\verb{[integer](1L) | [character](1L)}]\cr
Symbol used for indicating the local optima of the first objective
(default: \code{21L}, i.e., a filled circle).}

\item{symbol2}{[\verb{[integer](1L) | [character](1L)}]\cr
Symbol used for indicating the local optima of the second objective
(default: \code{22L}, i.e., a filled square).}

\item{symbol3}{[\verb{[integer](1L) | [character](1L)}]\cr
Symbol used for indicating the local optima of the third objective
(default: \code{24L}, i.e., a filled triangle).}

\item{...}{\link{any}\cr
Further arguments to be passed to the \code{geom_tile} function of \code{ggplot}.}
}
\value{
\code{\link{ggplot}}\cr
A \code{ggplot} object displaying the multi-objective gradient landscape.
}
\description{
If the underlying single-objective functions are \code{smoof}-functions,
extract their local optima and add them to the gradient field heatmap.
}
